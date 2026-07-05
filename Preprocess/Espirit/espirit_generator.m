function [sense_maps] = espirit_generator(kspace,ksize,ncalib,thr1,thr2)
%% ESPIRiT Maps
% Derives the ESPIRiT operator.
% Arguments:
% kspace: Multi channel k-space data. Expected dimensions are (sx, sy, nc), where (sx, sy) are volumetric dimensions and (nc) is the channel dimension.
% ksize: Parameter that determines the k-space kernel size. If kspace has dimensions (256, 256, 8), then the kernel will have dimensions (1, ksize, ksize, 8)
% ncalib: Parameter that determines the calibration region size. If kspace has dimensions (256, 256, 8), then the calibration region will have dimensions (ncalib, ncalib, 8)
% thr1: Parameter that determines the rank of the auto-calibration matrix (A). Singular values below thr1 times the largest singular value are set to zero.
% thr2: Crop threshold that determines eigenvalues "=1".
% Returns:
% sense_maps: This is the ESPIRiT operator. It will have dimensions (sx, sy, nc, nc) with (sx, sy, :, idx) being the idx'th set of ESPIRiT maps.

% get the sizes of the kspace data
[sx,sy,Nc] = size(kspace);

% crop a calibration area
calib = crop(kspace,[ncalib,ncalib,Nc]);

% Compute ESPIRiT EigenVectors
[k,S] = dat2Kernel(calib,[ksize ksize]);
idx = max(find(S >= S(1)*thr1));

% crop kernels and compute eigen-value decomposition in image space to get the maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

% crop sensitivity maps 
sense_maps = M(:,:,:,end).*repmat(W(:,:,end)>thr2,[1,1,Nc]);