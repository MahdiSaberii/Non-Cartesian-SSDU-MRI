load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')


load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_sbj6_2D.mat')
traj1 = trj_kxy;
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sbj6_2D.mat')
traj2 = trj_kxy;

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')


%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF_full = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
    par.FT{ie} = NUFFTDCF(traj1(:,:,1,ie), ones(size(DCF_full)), shift, [120, 120]);
end

shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')
for ie = 1:Neco
    DCF_10x = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{ie+6} = NUFFTDCF(traj2(:,1,ie), ones(size(DCF_10x)), shift, [120, 120]);
end

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

