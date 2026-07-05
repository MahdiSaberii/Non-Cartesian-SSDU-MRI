load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')
load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

%% %
load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcSingle_3_3_R8E12.mat')
dc10x = permute(dcf, [2,1]);

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcFull_3_3_R8E12.mat')
dcFull = permute(dcf, [3,2,1]);
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R8Echo12/full_scan2_sb6.mat')

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 12;
R = 8;
fprintf('Calculating NUFFT operators \n')
num = 1488;

for ie = 1:Neco
    DCF = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,:,1,ie), ones(size(DCF.^2)), shift, [80, 80]);
end
Ncoil = size(data_full,4); %different
Npar = 24;
img_uncombined_full = zeros([[80,80], Ncoil, Neco, Npar]);
data_full_reshaped = reshape(data_full, [num, Neco, R, Npar, Ncoil]);

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
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3/R8Echo12/single_scan2_sb6.mat')

%% Calculate NUFFT operator
shift   = [0,0];
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,1,ie), ones(size(DCF.^2)), shift, [80, 80]);
end
% 
Ncoil = size(data_single, 3); %different
Npar = 24;
img_uncombined = zeros([[80,80], Ncoil, Neco, Npar]);
data_single_reshaped = reshape(data_single, [num, Neco, 24, Ncoil]);

%%
fprintf('Gridding of k-space data \n')
sb = 0;%
count = 0;%%
for ip = 1:Npar
    coil_2 = permute(single(squeeze(coil(:,:,ip,:))), [3,1,2]);
    atb_2 = zeros(80,80,Neco);
    ref_2 = zeros(80,80,Neco);
    for ie = 1:Neco
        DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_single_reshaped(:,ie,ip,ic)), par.prec));  
        end
        atb_2(:,:,ie) = sum(permute(img_uncombined(:,:,:,ie,ip), [3,1,2]).*conj(coil_2), 1);
        ref_2(:,:,ie) = sum(permute(img_uncombined_full(:,:,:,ie,ip), [3,1,2]).*conj(coil_2), 1);
    end
    
    count = count + 1;
    
    atb_2 = permute(atb_2, [3,1,2]);
    ref_2 = permute(ref_2, [3,1,2]);
    idx = count;
    
    fname = sprintf('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_2/test/slice_%d.mat', count);
    save(fname, 'atb_2', 'ref_2', 'coil_2', 'idx', 'sb')
end
