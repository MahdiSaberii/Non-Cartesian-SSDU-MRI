import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import time


def fft_np(x):
    return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x, axes=[0,1]),axes=[0,1] , norm='ortho'), axes=[0,1]) 
def ifft_np(x):
    return np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(X , axes=[0,1]),axes=[0,1] ,norm='ortho'), axes=[0,1])

def espirit( ksp , ks = 6 , cs = 24 , t = 0.01 , c = 0.95 ):  # (kspace, kernel_size = 6 , calibration_size = 24 , t = 0.01 , c = 0.95)
    
    # Arguments 
    # cs= calibration size ,, t= determines the rank of matrix A ,, c= threshold for eigenvalues to be zero
    sy = ksp.shape[0]
    sz = ksp.shape[1]
    nc = ksp.shape[2]
    
    sy_calib = (sy//2-cs//2, sy//2+cs//2) if (sy > 1) else (0, 1)
    sz_calib = (sz//2-cs//2, sz//2+cs//2) if (sz > 1) else (0, 1)
    cksp = ksp[ sy_calib[0] : sy_calib[1] , sz_calib[0] : sz_calib[1] , : ]    #Cropped kspace    
    n_row = cksp.shape[0]
    n_column = cksp.shape[1]

    A = np.zeros([ (cs - ks +1)**2 , ( ks**2 ) * nc ] , dtype = complex)   # 2 is the number of dimension for kspace
    index = 0
    for row in range(n_row):
        for column in range(n_column):
            if ( (row + ks <= n_row) and (column + ks <= n_column) ):
                temp = cksp[ row : row + ks , column : column + ks , : ]
                temp = temp.flatten()
                A[ index , : ] = temp
                index = index + 1
    
    
    U, S, VH = np.linalg.svd(A, full_matrices=True)
    V = VH.conj().T
    n = np.sum(S >= t * S[0])
    V = V[:, 0:n]
    
    kernels = np.zeros(np.append(np.shape(ksp), n)).astype(np.complex64)
    kerdims = [ ks , ks , nc]
    
    kyt = (sy//2-ks//2 + 1, sy//2+ks//2 + 1) if (sy > 1) else (0, 1)
    kzt = (sz//2-ks//2 + 1, sz//2+ks//2 + 1) if (sz > 1) else (0, 1)
    for idx in range(n):
        kernels[ kyt[0] : kyt[1] , kzt[0] : kzt[1] , :, idx] = np.reshape(V[:, idx], kerdims)
    
    
    kerimgs = np.zeros(np.append(np.shape(ksp), n)).astype(np.complex64)
    for idx in range(n):
        for jdx in range(nc):
            ker = kernels[ ::-1 , ::-1 , jdx, idx].conj()
            kerimgs[:,:,jdx,idx] = fft_np(ker) * np.sqrt(sy * sz)/np.sqrt(ks**2)
    
    # Take the point-wise eigenvalue decomposition and keep eigenvalues greater than c
    maps = np.zeros(np.append(np.shape(ksp), nc)).astype(np.complex64)
    for jdx in range(sy):
        for kdx in range(sz):
            Gq = kerimgs[jdx,kdx,:,:]
            u, s, vh = np.linalg.svd(Gq, full_matrices=True)
            for ldx in range(nc):
                if (s[ldx]**2 > c):
                    maps[ jdx, kdx, :, ldx] = u[:, ldx]
    return maps # As I remember, the last dim is not necessary. So I took maps[:,:,:,0] for example for slices with the shape of [320,368,15]. Visualizing the maps helps a lot to take which dimension 


# kspace = sio.loadmat("./slice_66.mat")["kspace"] 
# maps = espirit(kspace)
# plt.imshow(np.abs(maps[:,:,3,0]), cmap="gray") # Last index should be zero

if __name__ == '__main__':
    
    DATASET = "PDFS_300"
    DATA_PATH = f"/home/naxos2-raid25/saber032/Main_works/Dataset/{DATASET}/Cropped/Test/Kspaces"
    SAVING_PATH = f"/home/naxos2-raid25/saber032/Main_works/Dataset/{DATASET}/Cropped/Test/Coils"
    print(f"{DATASET}")
    for i in range(0,395):
        ksp = sio.loadmat(f"{DATA_PATH}/slice_{i+1}.mat")["kspace"]
        tic = time.time()
        maps = espirit(ksp, t=0.003)[:,:,:,0]
        toc = time.time()
        print("Coil Number: ", i+1, "   Passed Time:", toc-tic)
        sio.savemat(f"{SAVING_PATH}/coil_{i+1}.mat", {"maps":maps})