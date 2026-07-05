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
traj1 = trj_kxy;

%%
%load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_sb6.mat')
%traj2 = trj_kxy;

%%
% a = traj1(:,:,1,1);
% b = traj2(:,:,1,1);

%%
% figure()
% plot(abs(a(:)))
% 
% figure()
% plot(abs(b(:)))
% 
% figure()
% plot(angle(a(:)))
% 
% figure()
% plot(angle(b(:)))

%%
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sbj6_2D.mat')
% traj1 = trj_kxy;
% 
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb6.mat')
% traj2 = trj_kxy;
% 
% %%
% a = traj1(:,1,1);
% b = traj2(:,1,1);
% 
% %%
% figure()
% plot(abs(a(:)))
% 
% figure()
% plot(abs(b(:)))
% 
% figure()
% plot(angle(a(:)))
% 
% figure()
% plot(angle(b(:)))

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
%Neco = par.Neco;
%Ncoil = par.NcoilFinal;
%Npar = par.Npar;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    DCF = dcFull(:,:,ie)./ max(max(dcFull(:,:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
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
Ncoil = size(data_single,4); %different
Npar = 24;
img_uncombined_full = zeros([[120,120], Ncoil, 6, Npar]);
data_full_reshaped = reshape(data_single, [2176, 6, 10, Npar, Ncoil]);

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

%%
figure()
for i = 1:6
    subplot(3,6,i)
    imshow(squeeze(abs(ref_2(i,:,:))), [])
    
    subplot(3,6,i+6)
    imshow(squeeze(abs(atb_2(i,:,:))), [])
    
    subplot(3,6,i+12)
    imshow(squeeze(abs(recon1(i,:,:))), [])
end

%%
ref_2 = zeros(120,120,6);
coil_2 = squeeze(coil(:,:,6,:));
for ie = 1:6
    ref_2(:,:,ie) = sum(img_uncombined_full(:,:,:,ie,6).*conj(coil_2), 3);
end

%%
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/1_2_3_4_5_6/large/reconst5.mat'

%%
rec_1 = squeeze(reconstruct(31,:,:,:));
rec_50 = squeeze(reconstruct(81,:,:,:));
rec_175 = squeeze(reconstruct(205,:,:,:));

%%
figure()
for i = 1:6
    subplot(4,6,i)
    imshow(squeeze(abs(ref_2(:,:,i))), [])
    
    subplot(4,6,i+6)
    imshow(squeeze(abs(rec_1(i,:,:))), [])
    
    subplot(4,6,i+12)
    imshow(squeeze(abs(rec_50(i,:,:))), [])
    
    subplot(4,6,i+18)
    imshow(squeeze(abs(rec_175(i,:,:))), [])
end


%%
img = permute(img_uncombined_full, [1,2,5,3,4]);
%%
coil = ZcGetSpiralCoils_2D(img,1);
coil = coil(:,:,:,:,1);

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