load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(abs(dc10x(:,ie))));
    FT{ie} = NUFFTDCF(trj_kxy(:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
end

Ncoil = size(data_single, 3); %different
Npar = 72;
img_uncombined = zeros([[120,120], Ncoil, Neco, Npar]);
data_single_reshaped = reshape(data_single, [2176, 6, 72, Ncoil]);

%%
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:Neco
        DCF = dc10x(:,ie)./ max(max(abs(dc10x(:,ie))));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined(:,:,ic,ie,ip) = FT{ie}' * (DCF.*cast(squeeze(data_single_reshaped(:,ie,ip,ic)), par.prec));          
        end
    end
end

% save sb7 img_uncombined

