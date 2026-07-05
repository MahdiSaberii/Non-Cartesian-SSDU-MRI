import torch
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import time

def fft_torch(x):
    return torch.fft.fftshift(torch.fft.fftn(torch.fft.ifftshift(x, dim=[0,1]), dim=[0,1] , norm='ortho'), dim=[0,1]) 

def ifft_torch(x):
    return torch.fft.ifftshift(torch.fft.ifftn(torch.fft.fftshift(x, dim=[0,1]), dim=[0,1] , norm='ortho'), dim=[0,1])

def espirit(ksp, ks=6, cs=24, t=0.01, c=0.95):
    sy = ksp.shape[0]
    sz = ksp.shape[1]
    nc = ksp.shape[2]

    sy_calib = (sy//2-cs//2, sy//2+cs//2) if (sy > 1) else (0, 1)
    sz_calib = (sz//2-cs//2, sz//2+cs//2) if (sz > 1) else (0, 1)
    cksp = ksp[slice(*sy_calib), slice(*sz_calib), :]    
    n_row = cksp.shape[0]
    n_column = cksp.shape[1]

    A = torch.zeros(((cs - ks + 1) ** 2, ks ** 2 * nc), dtype=torch.complex64)

    index = 0
    for row in range(n_row):
        for column in range(n_column):
            if (row + ks <= n_row) and (column + ks <= n_column):
                temp = cksp[row: row + ks, column: column + ks, :]
                temp = temp.flatten()
                A[index, :] = temp
                index = index + 1

    U, S, VH = torch.linalg.svd(A)
    V = VH.conj().T
    n = torch.sum(S >= t * S[0])
    V = V[:, :n]

    kernels = torch.zeros((*ksp.shape, n), dtype=torch.complex64)
    kerdims = [ks, ks, nc]

    kyt = (sy//2-ks//2 + 1, sy//2+ks//2 + 1) if (sy > 1) else (0, 1)
    kzt = (sz//2-ks//2 + 1, sz//2+ks//2 + 1) if (sz > 1) else (0, 1)
    for idx in range(n):
        kernels[kyt[0]: kyt[1], kzt[0]: kzt[1], :, idx] = torch.reshape(V[:, idx], kerdims)

    kerimgs = torch.zeros((*ksp.shape, n), dtype=torch.complex64)
    for idx in range(n):
        for jdx in range(nc):
            ker = kernels[..., jdx, idx].conj()
            kerimgs[..., jdx, idx] = fft_torch(ker) * torch.sqrt(torch.tensor(sy * sz, dtype=torch.float32)) / torch.sqrt(torch.tensor(ks ** 2, dtype=torch.float32))

    maps = torch.zeros((*ksp.shape, nc), dtype=torch.complex64)
    for jdx in range(sy):
        for kdx in range(sz):
            Gq = kerimgs[jdx, kdx, ...]
            u, s, vh = torch.linalg.svd(Gq)
            for ldx in range(nc):
                if s[ldx] ** 2 > c:
                    maps[jdx, kdx, :, ldx] = u[:, ldx]

    return maps

if __name__ == '__main__':
    DATASET = "PDFS_300"
    DATA_PATH = f"/home/naxos2-raid25/saber032/Main_works/Dataset/{DATASET}/Cropped/Kspaces"
    SAVING_PATH = f"/home/naxos2-raid25/saber032/Main_works/Dataset/{DATASET}/Cropped/Coils"
    print(f"{DATASET}")
    device = torch.device("cuda:3")
    for i in range(100, 300):
        ksp = sio.loadmat(f"{DATA_PATH}/slice_{i+1}.mat")["kspace"]
        ksp_torch = torch.tensor(ksp, dtype=torch.complex64, device=device)
        tic = time.time()
        maps = espirit(ksp_torch, t=0.003)[..., 0]
        toc = time.time()
        print("Coil Number: ", i+1, "   Passed Time:", toc - tic)
        sio.savemat(f"{SAVING_PATH}/coil_{i+1}.mat", {"maps": maps.cpu().numpy()})
