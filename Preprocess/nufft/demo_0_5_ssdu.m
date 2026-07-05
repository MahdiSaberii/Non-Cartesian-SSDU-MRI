load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sbj6_2D.mat')

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/mask_new1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/mask_new2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/mask_new3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/mask_new4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/mask_new5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for m = 1: 5
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{m, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [120, 120]);
end
end

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
            img_uncombined_validation(m, :,:,ic,ie,ip) = par.FT{m, ie}' * (DCF.*cast(squeeze(validation_data(m, :,ie,ip,ic)), par.prec));
        end
    end
end
end

%%
% 
img_gridding_sos = squeeze(sqrt(sum(abs(squeeze(img_uncombined_train(2,:,:,:,:,:))).^2,3)));
im = permute(squeeze(img_gridding_sos),[1 2 4 3]);

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

