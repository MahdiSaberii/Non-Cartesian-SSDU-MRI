load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

%% %
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_cor_2D_sbj1.mat')

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
end
Ncoil = size(data_full,4); %different
Npar = 72;
img_uncombined_full = zeros([[120,120], Ncoil, Neco, Npar]);
data_full_reshaped = reshape(data_full, [2176, 6, 10, Npar, Ncoil]);

%%
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:Neco
        DCF = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined_full(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_full_reshaped(:,ie,:,ip,ic)), par.prec));          
        end
    end
end
img = permute(img_uncombined_full, [1,2,5,3,4]);
coil = ZcGetSpiralCoils_2D(img,1);
coil = coil(:,:,:,:,1);

%% %
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleshot_cor_2D_sbj1.mat')

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
end
% 
Ncoil = size(data_single, 3); %different
Npar = 72;
img_uncombined = zeros([[120,120], Ncoil, Neco, Npar]);
data_single_reshaped = reshape(data_single, [2176, 6, 72, Ncoil]);

%%
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:Neco
        DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_single_reshaped(:,ie,ip,ic)), par.prec));          
        end
    end
end


%%
fprintf('Old mask \n')
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
fprintf('Extra mask \n')
training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;

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

Ncoil = size(data_single, 3); %different
Npar = 72;
img_uncombined_train = zeros([Ncoil, [120,120],Neco]);
img_uncombined_validation = zeros([Ncoil, [120,120], Neco]);

training_data = zeros(5, 2176, 6, 72, Ncoil);
validation_data  = zeros(5,2176, 6, 72, Ncoil);

for m = 1: 5
    training_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    training_data(m, training_all(m,:)~=1,:,:,:) = 0;
    validation_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    validation_data(m, training_all(m,:)==1,:,:,:) = 0;
end

%%
atb_2 = zeros(120, 120, 6);
fprintf('Gridding of k-space data \n')
sb = 1;%
count = 2100;%%
for ip = 7:Npar-6
    coil_2 = permute(single(squeeze(coil(:,:,ip,:))), [3,1,2]);
    for m = 1: 5
    atb_2_train = zeros(120,120,6);
    atb_2_val = zeros(120,120,6);
    atb_2 = zeros(120, 120, 6);
    for ie = 1:Neco
        DCF = transpose(dc10x(:,ie)./ max(max(dc10x(:,ie))));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined_train(ic, :,:,ie) = par.FT{m, ie}' * (DCF.*cast(squeeze(training_data(m, :,ie,ip,ic)), par.prec));          
            img_uncombined_validation(ic, :,:,ie) = par.FT{m+5, ie}' * (DCF.*cast(squeeze(validation_data(m, :,ie,ip,ic)), par.prec));
        end
        atb_2_train(:,:,ie) = sum(img_uncombined_train(:,:,:,ie).*conj(coil_2), 1);
        atb_2_val(:,:,ie) = sum(img_uncombined_validation(:,:,:,ie).*conj(coil_2), 1);
        atb_2(:,:,ie) = sum(permute(img_uncombined(:,:,:,ie,ip), [3,1,2]).*conj(coil_2),1);
    end
    count = count + 1;
    
    atb_2 = permute(atb_2_train, [3,1,2]);
    refi_2 = permute(atb_2_val, [3,1,2]);
    idx = count;
    
    m = m + 5; %%
    fname = sprintf('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/slices_lambda/slice_%d.mat', count);
    save(fname, 'atb_2', 'refi_2', 'coil_2', 'idx', 'sb', 'm')
    
    ref_2 = permute(atb_2, [3,1,2]);
    
    fname = sprintf('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/slices_omega/slice_%d.mat', count);
    save(fname, 'atb_2', 'ref_2', 'coil_2', 'idx', 'sb', 'm')
    m = m - 5; %%
    
    end
end