import os
import torch
import random
import numpy as np
from DC import DC_2d
from tqdm import tqdm
import scipy.io as sio
import matplotlib.pyplot as plt
from Unrolled_Network import UnrolledNet
from torchvision.utils import save_image
from DataLoader_Test import DataLoaderSL as DL


seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False


if __name__ == "__main__":
    
    if torch.cuda.is_available():  
        dev = "cuda:3" 
    else:  
        dev = "cpu"  
    device = torch.device(dev)
    MODEL            = "ConvPDDL" # "ConvPDDL" or "UMPIRE" or "PELPF" or "PENN"
    batch_size       = 1
    
    R              = 6
    Number_Masks   = 3 # This determines which model to load
    nEcho          = 5
    img_width      = 120
    threshold      = 0.00001
    window_scale   = 0.6
    
#       Reading the Data with DataLoaderSL
    Images_path_theta = f"../Data/Ehy_{R}_fMRI/" 
    CartesianData = DL(Images_path_theta, nEcho=nEcho)
    Data_sampler  = torch.utils.data.RandomSampler(CartesianData)
    data_loader   = torch.utils.data.DataLoader(dataset=CartesianData,
                                            batch_size=batch_size, 
                                            sampler = Data_sampler,
                                            num_workers = 4)
    
    maskkk       = torch.tensor(sio.loadmat("../Data/Tuke_01.mat")["LPfilter"]).to(device)
    dataset_size = len(data_loader)
    

    # Learning Setup 
    if   MODEL == "ConvPDDL":
        network    = UnrolledNet(device = device, Unrolls=10, echoes=nEcho).to(device)
    elif MODEL == "PELPF":
        network    = UnrolledNet(device = device, Unrolls=10, echoes=nEcho).to(device)
    network.load_state_dict(torch.load(f"./{MODEL}_{Number_Masks}Masks/model/nework_r{R}_epoch95.pth"))
    network.train()

    Saving_Folder     = f"./{MODEL}_{Number_Masks}Masks_Test/"

    if not os.path.exists(Saving_Folder):
        os.makedirs(Saving_Folder, exist_ok=True)
        os.makedirs(f"{Saving_Folder}pngs", exist_ok=True)
        os.makedirs(f"{Saving_Folder}model", exist_ok=True)              
    
    cg         = DC_2d(mu=0, device=device, iters=15, echos=nEcho, flag_learnable=False)
    train_loss = []  
    mu_values  = [[] for _ in range(nEcho)]
    
    num_params = sum(p.numel() for p in network.parameters() if p.requires_grad)
    print(f"Trainable parameters: {num_params:,}")

    progress_bar = tqdm(data_loader, total=len(data_loader), desc=f"Test fMRI!")
    for idx_loop, (omega_EHWy, coils,  sliceNumber, subjectNumber, time_index, FileName, M_single, label_ehy, M_full) in enumerate(progress_bar): 
        
        coil = coils.permute(1,0,2,3).squeeze() # Original shape: (1,52,120,120) -> (52,120,120)
    
        # move everything to gpu
        omega_EHWy  = omega_EHWy.to(device)        # Shape: [1,5,120,120]
        coil        = coil.to(device)              # Shape: (52,120,120)
        label_ehy   = label_ehy.to(device)
        
        # Masking the outer circle, because of the corner artifacts (this could have done during the preprocessing)
        label_ehy     = label_ehy * maskkk
        used_M_full   = M_full.to(device)
        M_single      = M_single.to(device)
    
        used_M_theta    = M_single # M.shape = [1,N_echoes,img_width,img_width] --> [1,5,120,120]
        output          = network(omega_EHWy , coil , used_M_theta) # Shape: (1,nEcho,120,120)
    
        

        dc_out = cg(omega_EHWy , omega_EHWy ,coil, used_M_theta)

        c_sbj    = int(subjectNumber.item())
        s_number = int(sliceNumber.item())
        t_index  = int(time_index.item())

        output        [torch.abs(label_ehy)<threshold] = 0
        dc_out        [torch.abs(label_ehy)<threshold] = 0
        label_ehy    [torch.abs(label_ehy)<threshold] = 0
        omega_EHWy    [torch.abs(label_ehy)<threshold] = 0
        
        x = torch.cat([output[0, i, :, :].cpu()/ torch.max(torch.abs(output[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
        z = torch.cat([dc_out[0, i, :, :].cpu()/ torch.max(torch.abs(dc_out[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
        w = torch.cat([omega_EHWy[0, i, :, :].cpu()/torch.max(torch.abs(omega_EHWy[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
        y = torch.cat([label_ehy[0, i, :, :].cpu()/ torch.max(torch.abs(label_ehy[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
        
        xx = torch.cat((y,x,z,w), dim=0)
        xx = xx / torch.max(torch.abs(xx))
        xx = torch.abs(xx).cpu().squeeze().detach()
        xx = torch.clamp(xx, max=window_scale) / window_scale
    
        save_image(torch.abs(xx), f'{Saving_Folder}pngs/Slc_{s_number:02d}_Sbj_{c_sbj:02d}_Time_{t_index:03d}.png')

        progress_bar.set_postfix({
                f"Sbj_{subjectNumber.cpu().item()} Slc_{sliceNumber.cpu().item()},"
                "Iter": f"{idx_loop+1}/{dataset_size}"})
    
    
