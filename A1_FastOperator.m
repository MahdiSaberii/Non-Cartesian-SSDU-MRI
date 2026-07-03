close all; clear all; clc;
% addpath("./nufft/")

N_Masks             = 7;
nEcho               = 5;
nArms               = 6;
TotalNumberofPoints = 3592;
img_width           = 120;
shift               = [0,0];
save_path           = './M/';

if ~exist(save_path, 'dir')
    mkdir(save_path);
end
%%
training_all        = zeros(N_Masks,TotalNumberofPoints);
validation_all      = zeros(N_Masks, TotalNumberofPoints);

for i=1:N_Masks
    name = "./Masks/mask_"+num2str(i)+".mat";
    load(name)
    training_all(i, :)  = training;
    validation_all(i,:) = validation;
end

%% Calculate NUFFT operator
load('./DCF/DCF_Full.mat')
trj_kxy      = kx + 1j*ky;
mask_all_one = ifftshift(ifftshift(load("./all_ones_mask.mat").mask_M, 1),2); 

%% Supervised M full
working_dcf = reshape(dcf_full, [TotalNumberofPoints*nArms, nEcho]);
kxx         = reshape(kx, [TotalNumberofPoints*nArms, nEcho]);
kyy         = reshape(ky, [TotalNumberofPoints*nArms, nEcho]);
M           = zeros(nEcho, img_width*2, img_width*2);
trj_full    = kxx + 1j*kyy;

for ie = 1:nEcho
    trj       = trj_full(:, ie);
    DCF       = working_dcf(:,ie)./ max(max(working_dcf(:,ie)));
    FTT       = NUFFTDCF(trj,ones(size(DCF)) , shift, [img_width*2, img_width*2]);

    W1        = FTT' * DCF;
    M1        = fft(fft(fftshift(fftshift(W1,1),2),[],1),[],2)/sqrt(img_width*2*img_width*2);
    M(ie,:,:) = M1 + mask_all_one .* mean(abs(M1(:)));
end
save(sprintf("%sM_full.mat",save_path), "M")

%% Calculate NUFFT operator
load('./DCF/DCF_Single.mat')
trj_kxy = kx + 1j*ky;
%% SSDU
M     = zeros(N_Masks*2, nEcho, img_width*2, img_width*2);
for m = 1: N_Masks
    working_dcf = dcf_single;
    trj_train   = trj_kxy;
    trj_val     = trj_kxy;

    trj_train(training_all(m,:)~=1, :) = 0;
    trj_val(validation_all(m,:)~=1, :) = 0;

    for ie = 1:nEcho
        
        DCF_train = working_dcf(:,ie)./ max(max(working_dcf(:,ie)));
        DCF_val   = working_dcf(:,ie)./ max(max(working_dcf(:,ie)));
        
        DCF_train(training_all(m,:)~=1) = 0;
        DCF_val(validation_all(m,:)~=1) = 0;

        FT_train = NUFFTDCF(trj_train(:,ie),ones(size(DCF_train)) , shift, [img_width*2, img_width*2]);
        FT_val   = NUFFTDCF(trj_val(:,ie)  ,ones(size(DCF_val))   , shift, [img_width*2, img_width*2]);
        grid_train = FT_train' * DCF_train;
        grid_val   = FT_val'   * DCF_val;
        
        M_theta             = fft(fft(fftshift(fftshift(grid_train,1),2),[],1),[],2)/sqrt(img_width*2*img_width*2);
        M(m,ie,:,:)         = M_theta + mask_all_one .* mean(abs(M_theta(:)));
        M_lambda            = fft(fft(fftshift(fftshift(grid_val,1),2),[],1),[],2)/sqrt(img_width*2*img_width*2);
        M(m+N_Masks,ie,:,:) = M_lambda + mask_all_one .* mean(abs(M_lambda(:)));
        disp(['Generating M_Theta, Echo ', num2str(ie),'/5' , '  R=' , num2str(6), ' Mask: ', num2str(m)])
    end
end
save(sprintf("%sM_single.mat", save_path)  , "M");

%% M Test Generation
load('./DCF/DCF_Single.mat')
trj_kxy = kx + 1j*ky;
mask_all_one = ifftshift(ifftshift(load("./all_ones_mask.mat").mask_M, 1),2); 
M            = zeros(nEcho, img_width*2, img_width*2);
    
working_dcf = dcf_single;
trj_train   = trj_kxy;

for ie = 1:nEcho
        
    DCF_train  = working_dcf(:,ie)./ max(max(working_dcf(:,ie)));
    FT_train   = NUFFTDCF(trj_train(:,ie),ones(size(DCF_train)) , shift, [img_width*2, img_width*2]);
    grid_train = FT_train' * DCF_train;
        
    M_theta   = fft(fft(fftshift(fftshift(grid_train,1),2),[],1),[],2)/sqrt(img_width*2*img_width*2);
    M(ie,:,:) = M_theta + mask_all_one .* mean(abs(M_theta(:)));

    disp(['Generating M_Theta, Echo ', num2str(ie),'/5' , '  R=' , num2str(6)])
end
save(sprintf("%sM_Test.mat",save_path)  , "M");