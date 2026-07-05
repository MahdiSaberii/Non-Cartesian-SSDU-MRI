%% Demonstration of data structure
% Golden angle multiecho 3D RAVE data acquisition (GA 3D stack-of-stars) for fat water separation, literature:
% https://doi.org/10.1002/mrm.26392,  https://doi.org/10.1002/mrm.28280

% free breathing dataset, no data sorting/binning performed, Golden Angle
% rotation between radial spokes => retrospective undersampling possible by
% reducing line dimension eg. reduced dataset rkdata = mkdata(:,1:end/2,:,:,:)
% modifications performed on kdata: FT along partitions, coilcompression, sampling density compensation
% modified kdata:            mkdata (#points per radial line, #radial lines, #coils, #echoes, #z partitions)
% trajectory coordinates:    traj (#points per radial line,#radial lines, #echoes, #z partitions)
% density compenstation function: dcf

load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')
training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;
%%
DCF_training_all = zeros(5, 1306, 6);
DCF_validation_all = zeros(5, 870, 6);
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_1.mat')
DCF_training_all(1, :, :) = dc10x_training;
DCF_validation_all(1, :, :) = dc10x_validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_2.mat')
DCF_training_all(2, :, :) = dc10x_training;
DCF_validation_all(2, :, :) = dc10x_validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_3.mat')
DCF_training_all(3, :, :) = dc10x_training;
DCF_validation_all(3, :, :) = dc10x_validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_4.mat')
DCF_training_all(4, :, :) = dc10x_training;
DCF_validation_all(4, :, :) = dc10x_validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_5.mat')
DCF_training_all(5, :, :) = dc10x_training;
DCF_validation_all(5, :, :) = dc10x_validation;
%%
%Neco = par.Neco;
%Ncoil = par.NcoilFinal;
%Npar = par.Npar;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating N operators \n')
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')



trj_training = zeros(5, 1306, 6);
trj_validation = zeros(5, 870, 6);

for m = 1: 5
for ie = 1:Neco
    trj_training(m,:,ie) = trj_kxy(training_all(m, :)==1, 1, ie);
    DCF = DCF_training_all(m, :,ie)./ max(max(DCF_training_all(m, :,ie)));
    par.FT{m, ie} = NUFFTDCF(trj_training(m,:,ie), DCF.^2, shift, [240, 240]);
end
end


for m = 1: 5
for ie = 1:Neco
    trj_validation(m,:,ie) = trj_kxy(validation_all(m, :)==1, 1, ie);
    DCF = DCF_validation_all(m, :,ie)./ max(max(DCF_validation_all(m, :,ie)));
    par.FT{m+5, ie} = NUFFTDCF(trj_validation(m,:,ie), DCF.^2, shift, [240, 240]);
end
end

% W_uncombined = zeros([par.imsize*4, Ncoil]);
% M_uncombined = zeros([par.imsize*4, Ncoil]);
% for ip = 1: Npar
%     W_uncombined(:,:,ip) = par.FT{1}' * dcf(:,:,ip);
%     M_uncombined(:,:,ip) = fftshift(fft(ifftshift(W_uncombined(:,:,ip)))) ./ (sqrt((512*2)^2));
% end

M = zeros(5, 6, 240, 240);
for m = 1: 5
for ie = 1:Neco
    W1 = par.FT{m, ie}' * (DCF_training_all(m,:,ie)/max(max(DCF_training_all(m,:,ie))));
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(m, ie,:,:) = M1;
end
end

%%
figure()
for i = 1:5
    temp = squeeze(M(i,1,:,:));
    subplot(2,3,i)
    imshow(abs(temp), [])

end
% W2 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1)))).^0.5;
% M2 = fftshift(fft2(ifftshift(W2))) ./ 240;
% 
% W3 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1))));
% M3 = fftshift(fft2(ifftshift(W3))) ./ 240;

% W1 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1)))).^0;
% M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
%%
% figure()
% imshow(abs(W1), [])
for m = 1: 5
figure()
for ie = 1: 6
    subplot(2,3,ie)
    imshow(squeeze(abs(M(m, ie,:,:))), [])
end
end
% figure()
% imshow(abs(M1), [])
% 
% figure()
% imshow(abs(W2), [])
% 
% figure()
% imshow(abs(M2), [])
% 
% figure()
% imshow(abs(W3), [])
% 
% figure()
% imshow(abs(M3), [])

% %% Gridding reconstruction
% under_mkdata = zeros(512, 600, 4, 6, 44);
% under_mkdata(:, 1:4:end, :, :, :) = mkdata(:, 1:4:end, :, :, :);
% 
% img_uncombined = zeros([par.imsize, Ncoil, Neco, Npar]);
% fprintf('Gridding of k-space data \n')
% for ip = 1:Npar
%     for ie = 1:Neco
%         for ic = 1:Ncoil
%             if par.verbose
%                 fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)
%             end            
%             img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * cast(squeeze(under_mkdata(:,:,ic,ie,ip)), par.prec);          
%         end
%     end
% end
% % img = img_uncombined(:,:,1,1,1);
% % 
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]);
% 
% 
% %%
% under_mkdata_1 = squeeze(under_mkdata(:,:,:,1,:));
% img_uncombined_1 = squeeze(img_uncombined(:,:,:,1,:));



%%
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]


%%
load /home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb7_coil.mat



