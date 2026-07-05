% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R6Echo5/full_scan6_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R6Echo5/full_scan12_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R6Echo5/single_scan6_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R6Echo5/single_scan12_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R8Echo7/full_scan5_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R8Echo7/single_scan5_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R10Echo8/full_scan4_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2/R10Echo8/single_scan4_sb2.mat')
% 
% 
% %%
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R6Echo10/full_scan3_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R6Echo10/single_scan3_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R8Echo12/full_scan2_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R8Echo12/single_scan2_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R10Echo15/full_scan1_sb2.mat')
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R10Echo15/single_scan1_sb2.mat')
% 
% %%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_full_R6E5.mat')
t1 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_single_R6E5.mat')
t2 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_full_R8E7.mat')
t3 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_single_R8E7.mat')
t4 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_full_R10E8.mat')
t5 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_single_R10E8.mat')
t6 = size(trj_kxy);

load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_full_R6E10.mat')
t7 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_single_R6E10.mat')
t8 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_full_R8E12.mat')
t9 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_single_R8E12.mat')
t10 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_full_R10E15.mat')
t11 = size(trj_kxy);
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_single_R10E15.mat')
t12 = size(trj_kxy);


%%
load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

%% %

for iter = 4:6
    if iter == 1
        file_string = '2_2/R6Echo5/';
        dcf_string = '2_2_R6E5';
        arm = t1(1);
        Necho = t1(4);
        R = t1(2);
        save_string = '2_2_1_1';
    elseif iter == 2
        file_string = '2_2/R8Echo7/';
        dcf_string = '2_2_R8E7';
        arm = t3(1);
        Necho = t3(4);
        R = t3(2);
        save_string = '2_2_2_1';
    elseif iter == 3
        file_string = '2_2/R10Echo8/';
        dcf_string = '2_2_R10E8';
        arm = t5(1);
        Necho = t5(4);
        R = t5(2);
        save_string = '2_2_3_1';
    elseif iter == 4
        file_string = '3_3/R6Echo10/'; 
        dcf_string = '3_3_R6E10';
        arm = t7(1);
        Necho = t7(4);
        R = t7(2);
        save_string = '3_3_1_1';
    elseif iter == 5
        file_string = '3_3/R8Echo12/';
        dcf_string = '3_3_R8E12';
        arm = t9(1);
        Necho = t9(4);
        R = t9(2);
        save_string = '3_3_2_1';
    else
        file_string = '3_3/R10Echo15/'; 
        dcf_string = '3_3_R10E15';
        arm = t11(1);
        Necho = t11(4);
        R = t11(2);
        save_string = '3_3_3_1';
    end
    
        
count = 0;%%   
for sbj = 1:4
    
    if sbj < 6 && iter == 1
        sb_num = sbj+1;
        scan_num = 6;
    end
    if sbj < 6 && iter == 2
        sb_num = sbj+1;
        scan_num = 5;
    end
    if sbj < 6 && iter == 3
        sb_num = sbj+1;
        scan_num = 4;
    end
    if sbj < 6 && iter == 4
        sb_num = sbj+1;
        scan_num = 3;
    end
    if sbj < 6 && iter == 5
        sb_num = sbj+1;
        scan_num = 2;
    end
    if sbj < 6 && iter == 6
        sb_num = sbj+1;
        scan_num = 1;
    end
    if sbj >= 6 && iter == 1
        sb_num = sbj-4;
        scan_num = 12;
    end
    if sbj >= 6 && iter == 2
        sb_num = sbj-4;
        scan_num = 11;
    end
    if sbj >= 6 && iter == 3
        sb_num = sbj-4;
        scan_num = 10;
    end
    if sbj >= 6 && iter == 4
        sb_num = sbj-4;
        scan_num = 9;
    end
    if sbj >= 6 && iter == 5
        sb_num = sbj-4;
        scan_num = 8;
    end
    if sbj >= 6 && iter == 6
        sb_num = sbj-4;
        scan_num = 7;
    end

full_load_string = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/' , file_string , 'full_scan', num2str(scan_num),'_sb' , num2str(sb_num) , '.mat'];
full_dcf_string = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcFull_' , dcf_string, '.mat'];

load(full_load_string)
load(full_dcf_string)

dcFull = permute(dcf, [3,2,1]);

%% Calculate NUFFT operator
shift   = [0,0];
Neco = Necho;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
end
Ncoil = size(data_full,4); %different
Npar = 24;

img_uncombined_full = zeros([[120,120], Ncoil, Neco, Npar]);
data_full_reshaped = reshape(data_full, [arm, Neco, R, Npar, Ncoil]);

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

%%
%coil_all = zeros(4,24,52,120,120);
%coil_all(1,:,:,:,:) = permute(coil(:,:,:,:), [3, 4, 1, 2]);

%% %
single_load_string = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/' , file_string , 'single_scan' , num2str(scan_num) , '_sb' , num2str(sb_num) , '.mat'];
single_dcf_string = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcSingle_' , dcf_string , '.mat'];
load(single_load_string)
load(single_dcf_string)

dc10x = permute(dcf, [2,1]);

%% Calculate NUFFT operator
shift   = [0,0];
Neco = Necho;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
end
% 
Ncoil = size(data_single, 3); %different
Npar = 24;
img_uncombined = zeros([[120,120], Ncoil, Neco, Npar]);
data_single_reshaped = reshape(data_single, [arm, Neco, Npar, Ncoil]);

%%
% fprintf('Gridding of k-space data \n')
% for ip = 1:Npar
%     for ie = 1:Neco
%         DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
%         for ic = 1:Ncoil
%             fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
%             img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_single_reshaped(:,ie,ip,ic)), par.prec));          
%         end
%     end
% end


%%
fprintf('mask \n')
training_all =zeros(10,arm);
validation_all = zeros(10, arm);

string1 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_1.mat'];
load(string1)
training_all(1, :) = training;
validation_all(1, :) = validation;

string2 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_2.mat'];
load(string2)
training_all(2, :) = training;
validation_all(2, :) = validation;

string3 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_3.mat'];
load(string3)
training_all(3, :) = training;
validation_all(3, :) = validation;

string4 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_4.mat'];
load(string4)
training_all(4, :) = training;
validation_all(4, :) = validation;

string5 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_5.mat'];
load(string5)
training_all(5, :) = training;
validation_all(5, :) = validation;

string6 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_6.mat'];
load(string6)
training_all(6, :) = training;
validation_all(6, :) = validation;

string7 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_7.mat'];
load(string7)
training_all(7, :) = training;
validation_all(7, :) = validation;

string8 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_8.mat'];
load(string8)
training_all(8, :) = training;
validation_all(8, :) = validation;

string9 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_9.mat'];
load(string9)
training_all(9, :) = training;
validation_all(9, :) = validation;

string10 = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_' , num2str(iter) , '_10.mat'];
load(string10)
training_all(10, :) = training;
validation_all(10, :) = validation;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = Necho;
fprintf('Calculating NUFFT operators \n')

%trj_training = zeros(5, 1306, 6);
%trj_validation = zeros(5, 870, 6);

for m = 1: 10
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{m, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [120, 120]);
end
end

for m = 1: 10
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{m+10, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [120, 120]);
end
end

Ncoil = size(data_single, 3); %different
Npar = 24;
img_uncombined_train = zeros([Ncoil, [120,120],Neco]);
img_uncombined_validation = zeros([Ncoil, [120,120], Neco]);

training_data = zeros(10, arm, Neco, Npar, Ncoil);
validation_data  = zeros(10,arm, Neco, Npar, Ncoil);

for m = 1: 10
    training_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    training_data(m, training_all(m,:)~=1,:,:,:) = 0;
    validation_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    validation_data(m, training_all(m,:)==1,:,:,:) = 0;
end

%%
fprintf('Gridding of k-space data \n')
sb = 0;%
for m = 1: 10
    for ip = 1:Npar
    coil_2 = permute(single(squeeze(coil(:,:,ip,:))), [3,1,2]);
    atb_2_train = zeros(120,120,Neco);
    atb_2_val = zeros(120,120,Neco);
%    atb_10 = zeros(120, 120, 6);
    for ie = 1:Neco
        DCF = transpose(dc10x(:,ie)./ max(max(dc10x(:,ie))));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined_train(ic, :,:,ie) = par.FT{m, ie}' * (DCF.*cast(squeeze(training_data(m, :,ie,ip,ic)), par.prec));          
            img_uncombined_validation(ic, :,:,ie) = par.FT{m+10, ie}' * (DCF.*cast(squeeze(validation_data(m, :,ie,ip,ic)), par.prec));
        end
        atb_2_train(:,:,ie) = sum(img_uncombined_train(:,:,:,ie).*conj(coil_2), 1);
        atb_2_val(:,:,ie) = sum(img_uncombined_validation(:,:,:,ie).*conj(coil_2), 1);
%        atb_10(:,:,ie) = sum(permute(img_uncombined(:,:,:,ie,ip), [3,1,2]).*conj(coil_2),1);
    end
    count = count + 1;
    
    atb_2 = permute(atb_2_train, [3,1,2]);
    refi_2 = permute(atb_2_val, [3,1,2]);
    idx = count;
    
    save_string_long = ['/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/',save_string,'/slice_',num2str(count),'.mat'];
    save(save_string_long, 'atb_2', 'refi_2', 'coil_2', 'idx', 'sb', 'm')
    
%     ref_2 = permute(atb_10, [3,1,2]);
%     
%     fname = sprintf('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/slices_omega_7/slice_%d.mat', count);
%     save(fname, 'atb_2', 'ref_2', 'coil_2', 'idx', 'sb', 'm')
    
    end
end

end
end



