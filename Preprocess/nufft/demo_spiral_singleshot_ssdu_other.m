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
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')

%%
%load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask5.mat')
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
fprintf('Calculating NUFFT operators \n')

trj_training = zeros(5, 1306, 6);
trj_validation = zeros(5, 870, 6);

% for m = 1: 5
% for ie = 1:Neco
%     trj_training(m,:,ie) = trj_kxy(training_all(m, :)==1, 1, ie);
%     DCF = DCF_training_all(m, :,ie);%./ max(max(DCF_training_all(m, :,ie)));
%     FT_1{m, ie} = NUFFTDCF(trj_training(m,:,ie), ones(size(DCF)), shift, [120, 120]);
% end
% end
% 
% for m = 1: 5
% for ie = 1:Neco
%     trj_validation(m,:,ie) = trj_kxy(validation_all(m, :)==1, 1, ie);
%     DCF = DCF_validation_all(m, :,ie);%./ max(max(DCF_validation_all(m, :,ie)));
%     FT_2{m, ie} = NUFFTDCF(trj_validation(m,:,ie), ones(size(DCF)), shift, [120, 120]);
% end
% end

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
Ncoil = size(data_single, 3); %different
Npar = 72;
img_uncombined_train = zeros([5, [120,120], Ncoil, Neco, Npar]);
img_uncombined_validation = zeros([5, [120,120], Ncoil, Neco, Npar]);
data_single_reshaped = reshape(data_single, [2176, 6, 72, Ncoil]);

training_data = zeros(5, 1306, 6, 72, Ncoil);
validation_data  = zeros(5, 870, 6, 72, Ncoil);

for m = 1: 5
    training_data(m, :,:,:,:) = data_single_reshaped(training_all(m, :)==1, :,:,:);
    validation_data(m, :,:,:,:) = data_single_reshaped(validation_all(m, :)==1, :,:,:);
end

%%
fprintf('Gridding of k-space data \n')
for m = 1: 5
for ip = 1:Npar
    for ie = 1:Neco
        trj_training(m,:,ie) = trj_kxy(training_all(m, :)==1, 1, ie);
        trj_validation(m,:,ie) = trj_kxy(validation_all(m, :)==1, 1, ie);
        DCF_train = DCF_training_all(m, :,ie)./ max(max(DCF_training_all(m, :,ie)));
        DCF_validation = DCF_validation_all(m, :,ie)./ max(max(DCF_validation_all(m, :,ie)));
        
        FT_1 = NUFFTDCF(trj_training(m,:,ie), ones(size(DCF_train)), shift, [120, 120]);
        FT_2 = NUFFTDCF(trj_validation(m,:,ie), ones(size(DCF_validation)), shift, [120, 120]);
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined_train(m, :,:,ic,ie,ip) = FT_1' * (DCF_train .* cast(squeeze(training_data(m, :,ie,ip,ic)), par.prec));          
            img_uncombined_validation(m, :,:,ic,ie,ip) = FT_2' * (DCF_validation .* cast(squeeze(validation_data(m, :,ie,ip,ic)), par.prec));
        end 
    end
end
end

%%
%img = img_uncombined(:,:,1,1,1);
%figure()
%imshow(abs(img), [])
%%
% 
img_gridding_sos_train = squeeze(sqrt(sum(abs(squeeze(img_uncombined_train(2,:,:,:,:,:))).^2,3)));
im_train = permute(squeeze(img_gridding_sos_train),[1 2 4 3]);

img_gridding_sos_validation = squeeze(sqrt(sum(abs(squeeze(img_uncombined_validation(2,:,:,:,:,:))).^2,3)));
im_validation = permute(squeeze(img_gridding_sos_validation),[1 2 4 3]);
% 
% 
% %%
% under_mkdata_1 = squeeze(under_mkdata(:,:,:,1,:));
% img_uncombined_1 = squeeze(img_uncombined(:,:,:,1,:));

%figure()
%imshow(abs(im(:,:,15,1)), [])

%%
a = im_train(:,:,34,1);
b = im_validation(:,:,34,1);

diff = a - b;

%%
c = squeeze(training_data(1,:,1,34,:));
d = squeeze(validation_data(1,:,1,34,:));

max_c = max(max(abs(c)));
max_d = max(max(abs(d)));

%%
max_a = max(max(abs(a)));
max_b = max(max(abs(b)));

max_diff = max(max(abs(diff)));

%%

figure()
imshow(abs(a), [0, 0.005])
figure()
imshow(abs(b), [0, 0.002])

%%
% figure()
% for ll = 1:6
%     echo_num = ll;
% for kk = 1:1
% %    title(['slice', num2str(34)])
%     subplot(2,3,ll)
%     
%     imshow(abs(im(:,:,34,echo_num)), [])
% %    title(['slice', num2str(kk*8)])
% end
% end
% 
% %%
% % img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% % im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]
% % 
% % load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/Training/Label/sb1/sb_1_slice_36.mat')
% % 
% % %%
% % idx = 6;
% % figure()
% % imshow(abs(im(:,:,36,idx)), [])
% % 
% % figure()
% % a = label_SingleSlice(:,:,:,idx);
% % b = complex(a(:,:,1), a(:,:,2));
% % imshow(abs(b), [])
% % 
% % %%
% % load ('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb1_coil.mat')
% % figure()
% % for ii = 1:58
% %     subplot(6,10,ii)
% %     imshow(angle(coil(:,:,36,ii)), [])  
% % end
