import glob
import os
import scipy.io as sio
import torch
from torch.utils.data import Dataset
from torch.utils.data import Dataset
import random


class DataLoaderSL(Dataset):    
    def __init__(self , datapath_theta, num_data = 2, nEcho = 5):

        self.nEcho = nEcho
        print(f" ------> Looking for Theat data in:  , {datapath_theta}")
        
        if num_data != 1:
            self.datapath_theta = glob.glob(os.path.join(datapath_theta , "*.mat"))
        else:
            all_files = glob.glob(os.path.join(datapath_theta, "*.mat"))
            specific_slice_files = [f for f in all_files if "subject_5_slice_18_" in f]
            # random.seed(42)  # You can use any integer value here
            self.datapath_theta = random.sample(specific_slice_files, num_data)
        print(" ------> number of Theta examples found: " , len(self.datapath_theta))
        
        
    def __getitem__(self , index):
      # Test Omega Data
      Omega_EHWy  = sio.loadmat(self.datapath_theta[index])[f"theta_ehy_6"].transpose([2,0,1])     # Original Data shape: [120,120,5]     == [img_x,img_y,echoes]
      
      FileName        = self.datapath_theta[index].split('/')[-1]    #full_cor_2D_sbj0_traj_1_slice_1.mat
      SliceNumber     = int(FileName.split('_')[-3])
      SubjectNumber   = int(FileName.split('_')[-5])
      TimeIndex       = int(FileName.split('_')[-1].split('.')[0])
      # Split takes: subject_1_slice_1_t_55 and gives: ['subject', '1', 'slice', '1', 't', '55']

      CoilAddress = f"../Data/Coils/{SubjectNumber}/slice_{SliceNumber}.mat"
      Coil        = sio.loadmat(CoilAddress)["sens_map"].transpose([2,0,1])   # Original Coil shape: [120,120,44] -> [44,120,120]
      
      
      # M Operator
      M_single = sio.loadmat(f"../Data/M/M_Test.mat")["M"] # [5,240,240] == [Echoes, 2*img_width, 2*img_width]
      M_single = torch.tensor(M_single)
      M_full   = sio.loadmat(f"../Data/M/M_full.mat")["M"]   # [5,240,240]    == [Echoes, 2*img_width, 2*img_width]
      M_full   = torch.tensor(M_full)
      
      label_path = f"../Data/Labels/"+ "subject_" + str(SubjectNumber) + "_slice_"+str(SliceNumber) + ".mat"
      label_ehy  = sio.loadmat(label_path)[f"label_ehy_temp"].transpose([2,0,1])  # Original Data shape: [120,120,5] 
      
      return Omega_EHWy, Coil , SliceNumber , SubjectNumber, TimeIndex, FileName, M_single,label_ehy, M_full
    
    def __len__(self):
        return len(self.datapath_theta)