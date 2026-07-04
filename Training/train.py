import os
import torch
import random
import numpy as np
from DC import DC_2d
from tqdm import tqdm
import scipy.io as sio
import matplotlib.pyplot as plt
from torchvision.utils import save_image
from DataLoader import DataLoaderSL as DL
from Unrolled_Network import UnrolledNet, UnrolledNet_PE

seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

def L1_L2_Loss(recon, label):
    loss = ( torch.norm(recon-label , p=2) / torch.norm( label, p=2) ) + ( torch.norm(recon-label , p=1) / torch.norm( label , p=1))
    return loss
		
def Cropping_center(input_tensor, crop_size = 40, crop_dims = [-2,-1]):
    fft_dcout_center = torch.zeros(input_tensor.shape, dtype=input_tensor.dtype, device=input_tensor.device)
    
    center_x = input_tensor.shape[crop_dims[0]] // 2
    center_y = input_tensor.shape[crop_dims[1]] // 2
    start_x  = center_x - crop_size // 2
    end_x    = center_x + crop_size // 2
    start_y  = center_y - crop_size // 2
    end_y    = center_y + crop_size // 2
    fft_dcout_center[..., start_x:end_x, start_y:end_y] = input_tensor[..., start_x:end_x, start_y:end_y]
    
    return fft_dcout_center

def M_Loss(image, coil , M):
    coil_images = image * coil
    size_x      = coil_images.shape[1] // 2
    size_y      = coil_images.shape[2] // 2
    
    a   = torch.nn.functional.pad(coil_images,[size_x,size_x,size_y,size_y])
    e   = torch.fft.fftn(torch.fft.fftshift(a, dim=[-2,-1]),dim=[-2,-1], norm='ortho')
    f   = e * M
    g   = torch.fft.ifftshift(torch.fft.ifftn(f ,dim=[-2,-1] ,norm='ortho'), dim=[-2,-1])
    h   = g[:,size_x:size_x*3,size_y:size_y*3]
    out = h
    # out = torch.sum(h*torch.conj(coil),axis=0) # Coil combined loss
    return out

if __name__ == "__main__": 
   
    device           = "cuda:2" if torch.cuda.is_available() else "cpu"
    MODEL            = "ConvPDDL" # "ConvPDDL" or "PELPF"
    batch_size       = 1
    learning_rate    = 5e-4
    scheduler_step   = 30
    scheduler_gamma  = 0.9
    epochs_number    = 100
    LPF_size         = 32
    
    R                = 6
    Number_Masks     = 7
    nEcho            = 5
    img_width        = 120
    threshold        = 0.00001
    window_scale     = 0.6
    
#       Reading the Data with DataLoaderSL
    Images_path_theta = f"../Data/Ehy_{R}/" 
    CartesianData = DL(Images_path_theta, nEcho=nEcho)
    Data_sampler  = torch.utils.data.RandomSampler(CartesianData)
    data_loader   = torch.utils.data.DataLoader(dataset=CartesianData,
                                            batch_size=batch_size, 
                                            sampler = Data_sampler,
                                            num_workers = 4)
    
    maskkk       = torch.tensor(sio.loadmat(f"../Data/Tuke_01.mat")["LPfilter"]).to(device)
    LPF          = torch.tensor(sio.loadmat(f"../Data/Tukey_{LPF_size}.mat")["LPfilter"]).to(device)
    dataset_size = len(data_loader)
    

    # Learning Setup 
    if   MODEL == "ConvPDDL":
        network    = UnrolledNet(device = device, Unrolls=10, echoes=nEcho).to(device)
    elif MODEL == "PELPF":
        network    = UnrolledNet_PE(device = device, Unrolls=10, echoes=nEcho).to(device)
    # network.load_state_dict(torch.load(f"./BestModel/nework_r6_epoch119.pth"))
    network.train()

    Saving_Folder     = f"./{MODEL}_{Number_Masks}Masks/"


    if not os.path.exists(Saving_Folder):
        os.makedirs(Saving_Folder, exist_ok=True)
        os.makedirs(f"{Saving_Folder}pngs", exist_ok=True)
        os.makedirs(f"{Saving_Folder}model", exist_ok=True)

    optimizer  = torch.optim.Adam(network.parameters(), lr = learning_rate)
    scheduler  = torch.optim.lr_scheduler.StepLR(optimizer, step_size=scheduler_step, gamma=scheduler_gamma)              
    
    CG         = DC_2d(mu=0, iters=15, echos=nEcho, flag_learnable=False).to(device)
    train_loss = []  
    mu_values  = [[] for _ in range(nEcho)]
    
    num_params = sum(p.numel() for p in network.parameters() if p.requires_grad)
    print(f"Trainable parameters: {num_params:,}")

    for epoch in range(1,epochs_number):
        losss = 0
        progress_bar = tqdm(data_loader, total=len(data_loader), desc=f"Epoch {epoch}/{epochs_number}")
        for idx_loop, (theta_EHWy  ,lambda_EHWy, coils,  sliceNumber, subjectNumber, time_index, FileName, M_single, label_ehy, M_full) in enumerate(progress_bar): 
            
            coil = coils.permute(1,0,2,3).squeeze()    # Original shape: (1,52,120,120) -> (52,120,120)
        
            # move everything to gpu
            theta_EHWy  = theta_EHWy.to(device)        # Shape: [1,5,120,120,1,5]
            lambda_EHWy = lambda_EHWy.to(device)       # Shape: [1,5,44,120,120,1,5]          
            coil        = coil.to(device)              # Shape: (52,120,120)
            label_ehy   = label_ehy.to(device)
            
            n_coils     = lambda_EHWy.shape[2]
            label_ehy   = label_ehy * maskkk
            used_M_full = M_full.to(device)
            M_single    = M_single.to(device)

            loss = 0.0
            for i_mask in range(Number_Masks):
                used_M_theta    = M_single[:,i_mask,:,:,:].squeeze(1) # M.shape = [1,N_Masks*2,N_echoes,120,120] --> [1,N_echoes,120,120]
                # The M is originally generated with 7 MAsks
                used_M_lambda   = M_single[:,i_mask+7,:,:,:].squeeze(1) # M.shape = [1,N_Masks*2,N_echoes,120,120] --> [1,N_echoes,120,120] 
                theta_EHWy_idx  = theta_EHWy[:,:,:,:,i_mask].squeeze(-1).squeeze(-2) 
                lambda_EHWy_idx = lambda_EHWy[:,:,:,:,:,i_mask].squeeze(-1).squeeze(-2) 

				dc_out           = CG(theta_EHWy , theta_EHWy  ,coil, used_M_theta) 
	            fft_dcout        = torch.fft.fftshift(torch.fft.fftn(torch.fft.fftshift(dc_out, dim=[-2,-1]), dim=[-2,-1], norm='ortho'), dim=[-2,-1])
	            fft_dcout_center = Cropping_center(fft_dcout, crop_size=LPF_size, crop_dims=[-2,-1]) * LPF
	            dc_out_filtered  = torch.fft.ifftshift(torch.fft.ifftn(torch.fft.ifftshift(fft_dcout_center, dim=[-2,-1]), dim=[-2,-1], norm='ortho'), dim=[-2,-1])
	            p_prime          = dc_out_filtered / (torch.sqrt(dc_out_filtered.real**2 + dc_out_filtered.imag**2) + 1e-12)
				
                if   MODEL == "ConvPDDL":
                    output          = network(theta_EHWy_idx , coil , used_M_theta)          # Shape: (1,nEcho,120,120)
                elif MODEL == "PELPF":
                    output          = network(theta_EHWy_idx , coil , used_M_theta, p_prime) # Shape: (1,nEcho,120,120)

                # Loss Calculation
                loss_image = torch.zeros([1, nEcho, n_coils, img_width, img_width], dtype = torch.complex64).to(device)
                for idx in range(nEcho):
                    loss_image[:,idx,:,:,:] = M_Loss(output[:,idx,:,:] , coil, used_M_lambda[:,idx,:,:])
                
                loss1 = L1_L2_Loss(loss_image , lambda_EHWy_idx) # Shapes: torch.Size([1, nEcho, nCoils, img_width, img_width])
                loss  = loss + loss1/Number_Masks
            
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()

            dc_out = CG(theta_EHWy_idx , theta_EHWy_idx ,coil, used_M_theta)

            cc       = [9,12,15,22]
            c_sbj    = int(subjectNumber.item())
            s_number = int(sliceNumber.item())
            t_index  = int(time_index.item())

            if s_number in cc:
                output        [torch.abs(label_ehy)<threshold] = 0
                dc_out        [torch.abs(label_ehy)<threshold] = 0
                label_ehy     [torch.abs(label_ehy)<threshold] = 0
                theta_EHWy_idx[torch.abs(label_ehy)<threshold] = 0
                
                x = torch.cat([output[0, i, :, :].cpu()/ torch.max(torch.abs(output[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
                z = torch.cat([dc_out[0, i, :, :].cpu()/ torch.max(torch.abs(dc_out[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
                w = torch.cat([theta_EHWy_idx[0, i, :, :].cpu()/torch.max(torch.abs(theta_EHWy_idx[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
                y = torch.cat([label_ehy[0, i, :, :].cpu()/ torch.max(torch.abs(label_ehy[0, i, :, :].cpu()))  for i in range(nEcho)],dim=1)
                
                xx = torch.cat((y,x,z,w), dim=0)
                xx = xx / torch.max(torch.abs(xx))
                xx = torch.abs(xx).cpu().squeeze().detach()
                xx = torch.clamp(xx, max=window_scale) / window_scale
            
                save_image(torch.abs(xx), f'{Saving_Folder}pngs/E{epoch}_Slc_{s_number:02d}_Sbj_{c_sbj:02d}_Time_{t_index:03d}.png')
            
            model_path = f"{Saving_Folder}model/nework_r{R}_epoch{epoch}.pth"   
            torch.save(network.state_dict(), model_path)

            losss = losss + (loss.item())
            average_loss = losss/(idx_loop+1)
            progress_bar.set_postfix({
                    f"Sbj_{subjectNumber.cpu().item()} Slc_{sliceNumber.cpu().item()},"
                    "Total Loss"  : f"{losss / dataset_size:.4f}",
                    "average_loss": f"{(losss/(idx_loop+1)):.4f}",
                    "Iter": f"{idx_loop+1}/{dataset_size}"})
        
        scheduler.step()
        for idx in range(nEcho):
            mu_values[idx].append(network.DataConsistency.mu[idx].cpu().detach().numpy())

        train_loss.append(losss/ dataset_size)
        muu1_numpy = np.array(mu_values)
        train_loss_mat = np.array(train_loss)

        sio.savemat(f'{Saving_Folder}/loss_R{R}.mat' , {"loss": train_loss_mat})
        sio.savemat(f"{Saving_Folder}/mu_R{R}.mat",   {"mu"  : muu1_numpy})

        # Plot and save loss
        plt.figure(figsize=(10, 6))
        plt.plot(train_loss)
        plt.title(f'Training Loss (R={R})')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.grid(True)
        plt.savefig(f'{Saving_Folder}/loss_R{R}.jpg', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot and save mu values      
        plt.figure(figsize=(10, 6))
        for idx in range(nEcho):
            plt.plot(np.array(mu_values[idx]), label=f'Echo {idx+1}')
        plt.title(f'Mu Values (R={R})')
        plt.xlabel('Epochs')
        plt.ylabel('Mu Value')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'{Saving_Folder}/mu_R{R}.jpg', dpi=300, bbox_inches='tight')
        plt.close()
