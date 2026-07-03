import numpy as np
import sigpy.mri
import scipy.io as sio

    # Parameters
TotalNumberofPoints = 3592
device_number       = 3
n_echos             = 5
n_arms              = 6
N_Masks             = 7
img_width           = 120
data_path           = f"./Raw/"
density_path        = f"./DCF/"
DATA_NAME_REF       = "trj_kxy_full"   
DATA_NAME           = "trj_kxy_single" 


# Full Trajectory Part
    # Traj size:       17960,6   (TotalNumberofPoints*n_echos, n_arms)
full_data_traj = sio.loadmat(data_path + f"{DATA_NAME_REF}.mat")["traj_full"] 
full_data_traj = full_data_traj.transpose([0,2,1])  # [TotalNumberofPoints, n_arms, n_echos]
kx_full        = np.transpose(np.real(full_data_traj) , [2,1,0])
ky_full        = np.transpose(np.imag(full_data_traj) , [2,1,0])
dcf_full       = np.zeros([n_echos , n_arms , TotalNumberofPoints] , dtype = np.float32)

for i in range(n_echos):
    new_traj        = np.zeros([n_arms , TotalNumberofPoints , 2] , dtype = np.float32)
    new_traj[:,:,0] = kx_full[i,:,:]
    new_traj[:,:,1] = ky_full[i,:,:]
    new_traj        = new_traj * (img_width*2)   # traj should be in the range of the image
    dcf_full[i,:,:] = sigpy.mri.pipe_menon_dcf(new_traj,device = device_number, max_iter=200,n=128, beta=8, width=4, show_pbar=True).get()
sio.savemat( f"{density_path}DCF_Full.mat" , {"dcf_full" : dcf_full.transpose([2,1,0]), "kx": np.real(full_data_traj) , "ky": np.imag(full_data_traj)})


# SingleShot i.e. R=6
    # Traj size:       17960,1 (TotalNumberofPoints*n_echos, 1)
single_data_traj = sio.loadmat(data_path + f"{DATA_NAME}.mat")["traj_single"]

kx_single  = np.transpose(np.real(single_data_traj) , [1,0])
ky_single  = np.transpose(np.imag(single_data_traj) , [1,0])
dcf_single = np.zeros([n_echos , TotalNumberofPoints] , dtype = np.float32)

for i in range(n_echos):
    new_traj        = np.zeros([TotalNumberofPoints , 2] , dtype = np.float32)
    new_traj[:,0]   = kx_single[i,:]
    new_traj[:,1]   = ky_single[i,:]
    new_traj        = new_traj * (img_width*2)   # traj should be in the range of the image
    dcf_single[i,:] = sigpy.mri.pipe_menon_dcf(new_traj,device = device_number, max_iter=200,n=128, beta=8, width=4, show_pbar=True).get()
sio.savemat( f"{density_path}DCF_Single.mat" , {"dcf_single" : dcf_single.transpose([1,0]) , "kx": np.real(single_data_traj) , "ky": np.imag(single_data_traj)})
