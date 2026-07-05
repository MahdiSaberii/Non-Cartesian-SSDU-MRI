
import scipy.io as sio
import glob
from torch.utils.data import Dataset
import os
import numpy as np
import sigpy.mri






ksp = sio.loadmat("slice_19.mat")["kspace"]
mask = sio.loadmat("Omega_Mask.mat")["mask"]
# maps = sio.loadmat("slice_19_coil.mat")["maps"]
ksp = ksp.transpose([2,0,1])
# maps = maps.transpose([2,0,1])
maps = sigpy.mri.app.EspiritCalib(ksp, calib_width=24, thresh=0.005).run()
sio.savemat("./slice_19_coil.mat", {"maps":maps.transpose([1,2,0])})
zf = np.sum(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ksp*mask, axes=[1,2]), norm="ortho", axes=[1,2]),axes=[1,2])* np.conj(maps),axis=0)
label = np.sum(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ksp, axes=[1,2]), norm="ortho", axes=[1,2]),axes=[1,2])* np.conj(maps),axis=0)
sio.savemat("./slice_19_zf.mat", {"zf":zf})
sio.savemat("./slice_19_label.mat", {"label":label})

