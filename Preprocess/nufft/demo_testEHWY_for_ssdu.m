load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sbj6_2D.mat')
traj = trj_kxy;
%%
traj_singleshot = squeeze(traj(:,1,:));

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/SSDU725/Density/DCF_new1.mat')
%load('/home/naxos2-raid7/hongygu/non_cartesian/SSDU725/Mask_generation/mask_new1.mat')
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/mask_new1.mat')


load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

dc10x_training = dc10x(training == 1,:);
dc10x_validation = dc10x(validation == 1,:);

%%
traj_train = traj(training == 1,1,:);
traj_validation = traj(validation == 1,1,:);

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    curr_traj = traj_train(:,:,ie);
    par.FT{ie} = NUFFTDCF(curr_traj, ones(size(curr_traj)), shift, [120, 120]);
end

shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')
for ie = 1:Neco
    curr_traj = traj_validation(:,:,ie);
    par.FT{ie+6} = NUFFTDCF(curr_traj, ones(size(curr_traj)), shift, [120, 120]);
end

%%
slice = reshape(squeeze(data_single(:, 32, :)), [2176, 6, 52]);

slice_train = slice(training == 1,1,1);
slice_validation = slice(validation == 1,1,1);

DCF_train = dc10x_training(:,1:2).';
DCF_train = DCF_train(1,:);
DCF_train = DCF_train.';

DCF_validation = dc10x_validation(:,1:2).';
DCF_validation = DCF_validation(1,:);
DCF_validation = DCF_validation.';

full_rec_wu = par.FT{1}' * (cast(slice_train, par.prec));
x_rec_wu = par.FT{7}' * (cast(slice_validation, par.prec));

a = [full_rec_wu, x_rec_wu];

fraction1 = mean(abs(full_rec_wu))/mean(abs(x_rec_wu));
fraction2 = mean(abs(DCF_train))/mean(abs(DCF_validation));

figure()
imshow(abs(a), [])

full_rec_wu = par.FT{1}' * (DCF_train.*cast(slice_train, par.prec));
x_rec_wu = par.FT{7}' * (DCF_validation.* cast(slice_validation, par.prec));

a = [full_rec_wu, x_rec_wu];

fraction3 = mean(abs(full_rec_wu))/mean(abs(x_rec_wu));
fraction4 = mean(abs(DCF_train))/mean(abs(DCF_validation));

figure()
imshow(abs(a), [])

DCF_train_1 = DCF_train/max(abs(DCF_train(:)));
DCF_validation_1 = DCF_validation/max(abs(DCF_validation(:)));

full_rec_wu = par.FT{1}' * (DCF_train_1.*cast(slice_train, par.prec));
x_rec_wu = par.FT{7}' * (DCF_validation_1.* cast(slice_validation, par.prec));

a = [full_rec_wu, x_rec_wu];

fraction5 = mean(abs(full_rec_wu))/mean(abs(x_rec_wu));
fraction6 = mean(abs(DCF_train_1))/mean(abs(DCF_validation_1));

figure()
imshow(abs(a), [])

%%
data_full_reshaped = reshape(data_full, [2176, 6, 10, 72, 52]);
full_rec_wo = par.FT{1}' * (cast(squeeze(data_full_reshaped(:,1,:,34,1)), par.prec));
x_rec_wo = par.FT{7}' * (cast(squeeze(data_single(1:2176,34,1)), par.prec));

a = [full_rec_wo, x_rec_wo];

fraction1 = mean(abs(full_rec_wo))/mean(abs(x_rec_wo));

figure()
imshow(abs(a), [])

full_rec_wu = par.FT{1}' * (dcFull(:,:,ie).* cast(squeeze(data_full_reshaped(:,1,:,34,1)), par.prec));
x_rec_wu = par.FT{7}' * (dc10x(:,ie).* cast(squeeze(data_single(1:2176,34,1)), par.prec));

a = [full_rec_wu, x_rec_wu];

fraction2 = mean(abs(full_rec_wu))/mean(abs(x_rec_wu));

figure()
imshow(abs(a), [])

full_rec_wn = par.FT{1}' * (DCF_full.* cast(squeeze(data_full_reshaped(:,1,:,34,1)), par.prec));
x_rec_wn = par.FT{7}' * (DCF_10x.* cast(squeeze(data_single(1:2176,34,1)), par.prec));

a = [full_rec_wn, x_rec_wn];

fraction3 = mean(abs(full_rec_wn))/mean(abs(x_rec_wn));

figure()
imshow(abs(a), [])
% W_uncombined = zeros([par.imsize*4, Ncoil]);
% M_uncombined = zeros([par.imsize*4, Ncoil]);
% for ip = 1: Npar
%     W_uncombined(:,:,ip) = par.FT{1}' * dcf(:,:,ip);
%     M_uncombined(:,:,ip) = fftshift(fft(ifftshift(W_uncombined(:,:,ip)))) ./ (sqrt((512*2)^2));
% end
%%
% W = par.FT{1}' * dcf(:,:,1)^0;
% M = fftshift(fft(ifftshift(W_uncombined(:,:,ip)))) ./ 1024;
% 
% figure()
% imshow(W, [])
% 
% figure()
% imshow(M, [])

% %% Gridding reconstruction
% under_mkdata = zeros(512, 600, 4, 6, 44);
% under_mkdata(:, 1:4:end, :, :, :) = mkdata(:, 1:4:end, :, :, :);
% 
Ncoil = size(data_full, 4); %different
Npar = size(data_full,3);
img_uncombined = zeros([[120,120], Ncoil, Neco, Npar]);
data_full_reshaped = reshape(data_full, [2176, 6, 10, Npar, Ncoil]);

%%
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:Neco
        DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_full_reshaped(:,ie,1,ip,ic)), par.prec));          
        end
    end
end

%%
%img = img_uncombined(:,:,1,1,1);
%figure()
%imshow(abs(img), [])
%%
% 
img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
im = permute(squeeze(img_gridding_sos),[1 2 4 3]);
% 
% 
% %%
% under_mkdata_1 = squeeze(under_mkdata(:,:,:,1,:));
% img_uncombined_1 = squeeze(img_uncombined(:,:,:,1,:));

%figure()
%imshow(abs(im(:,:,15,1)), [])

%%
for ll = 1:6
    figure()
    echo_num = ll;
for kk = 1:8
    title(['slice', num2str(kk*8-8)])
    subplot(2,4,kk)
    
    imshow(abs(im(:,:,kk*8,echo_num)), [])
    title(['slice', num2str(kk*8)])
end
end

%%
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]

