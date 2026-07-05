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
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sbj6_2D.mat')

%%
%load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/0_4/mask_new_0_4_1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/0_4/mask_new_0_4_2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/0_4/mask_new_0_4_3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/0_4/mask_new_0_4_4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/0_4/mask_new_0_4_5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%

% DCF_training_all = zeros(5, 1306, 6);
% DCF_validation_all = zeros(5, 870, 6);
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_new1.mat')
% DCF_training_all(1, :, :) = dc10x_training;
% DCF_validation_all(1, :, :) = dc10x_validation;
% 
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_new2.mat')
% DCF_training_all(2, :, :) = dc10x_training;
% DCF_validation_all(2, :, :) = dc10x_validation;
% 
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_new3.mat')
% DCF_training_all(3, :, :) = dc10x_training;
% DCF_validation_all(3, :, :) = dc10x_validation;
% 
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_new4.mat')
% DCF_training_all(4, :, :) = dc10x_training;
% DCF_validation_all(4, :, :) = dc10x_validation;
% 
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/DCF/DCF_new5.mat')
% DCF_training_all(5, :, :) = dc10x_training;
% DCF_validation_all(5, :, :) = dc10x_validation;


%%
%Neco = par.Neco;
%Ncoil = par.NcoilFinal;
%Npar = par.Npar;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

%trj_training = zeros(5, 1306, 6);
%trj_validation = zeros(5, 870, 6);

for m = 1: 5
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{m, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [120, 120]);
end
end

for m = 1: 5
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{m+5, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [120, 120]);
end
end

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

training_data = zeros(5, 2176, 6, 72, Ncoil);
validation_data  = zeros(5,2176, 6, 72, Ncoil);

for m = 1: 5
    training_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    training_data(m, training_all(m,:)~=1,:,:,:) = 0;
    validation_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    validation_data(m, training_all(m,:)==1,:,:,:) = 0;
end

%%
fprintf('Gridding of k-space data \n')
for m = 1: 5
for ip = 1:Npar
    for ie = 1:Neco
        DCF = transpose(dc10x(:,ie)./ max(max(dc10x(:,ie))));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined_train(m, :,:,ic,ie,ip) = par.FT{m, ie}' * (DCF.*cast(squeeze(training_data(m, :,ie,ip,ic)), par.prec));          
            img_uncombined_validation(m, :,:,ic,ie,ip) = par.FT{m+5, ie}' * (DCF.*cast(squeeze(validation_data(m, :,ie,ip,ic)), par.prec));
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
img_gridding_sos = squeeze(sqrt(sum(abs(squeeze(img_uncombined_train(2,:,:,:,:,:))).^2,3)));
im = permute(squeeze(img_gridding_sos),[1 2 4 3]);
% 
% 
% %%
% under_mkdata_1 = squeeze(under_mkdata(:,:,:,1,:));
% img_uncombined_1 = squeeze(img_uncombined(:,:,:,1,:));

%figure()
%imshow(abs(im(:,:,15,1)), [])

%%
figure()
for ll = 1:6
    echo_num = ll;
for kk = 1:1
%    title(['slice', num2str(34)])
    subplot(2,3,ll)
    
    imshow(abs(im(:,:,34,echo_num)), [])
%    title(['slice', num2str(kk*8)])
end
end

%%
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]
% 
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/Training/Label/sb1/sb_1_slice_36.mat')
% 
% %%
% idx = 6;
% figure()
% imshow(abs(im(:,:,36,idx)), [])
% 
% figure()
% a = label_SingleSlice(:,:,:,idx);
% b = complex(a(:,:,1), a(:,:,2));
% imshow(abs(b), [])
% 
% %%
% load ('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb1_coil.mat')
% figure()
% for ii = 1:58
%     subplot(6,10,ii)
%     imshow(angle(coil(:,:,36,ii)), [])  
% end