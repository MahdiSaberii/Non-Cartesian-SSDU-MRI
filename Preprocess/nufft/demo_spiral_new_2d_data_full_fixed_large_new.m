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

close all
load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')

%%
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/large_testing/fullsample_2D_test0803.mat'

%%
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/large_testing/full_1023_ISMRM.mat'
traj1 = trj_kxy;
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
end

% 
Ncoil = size(data_full,4); %different
Npar = 24;
img_uncombined_full = zeros([[120,120], Ncoil, 6, Npar]);
data_full_reshaped = reshape(data_full, [2176, 6, 10, Npar, Ncoil]);

%%
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:6
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
combined = zeros(6,120,120);
for ie = 1:6
    combined(ie,:,:) = sum(img_uncombined_full(:,:,:,ie,12).*conj(squeeze(coil(:,:,12,:))), 3);
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
% for ll = 1:6
%     figure()
%     echo_num = ll;
% for kk = 1:8
%     title(['slice', num2str(kk*8-8)])
%     subplot(2,4,kk)
%     
%     imshow(abs(im(:,:,kk*8,echo_num)), [])
%     title(['slice', num2str(kk*8)])
% end
% end
%%
figure()
imshow(abs(im(:,:,34,1)), [])
%%
coil_1 = permute(coil, [4,1,2,3]);
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/new_data_full_sb6.mat')

a = permute(squeeze(img_uncombined(:,:,:,1,34)), [3,1,2]);
b = squeeze(coil_1(:,:,:,34));

result = a.*conj(b);

final = squeeze(sum(result, 1));

figure()
imshow(abs(final), [])
% %%
% % img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% % im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]
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
%%
a = data_full;
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_sbj2_2D.mat')

b = a - data_full;