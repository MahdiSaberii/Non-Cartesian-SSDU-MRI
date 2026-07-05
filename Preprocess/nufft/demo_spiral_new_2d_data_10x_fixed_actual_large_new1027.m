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

%close all

load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleshot_cor_2D_sbj0.mat'
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/large_testing/fmri_1027_ISMRM.mat'
traj1 = trj_kxy;
%load '/home/naxos2-raid7/hongygu/non_cartesian/coil_large.mat'

%%
%load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/large_testing/fullsample_2D_test0803.mat'
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/coil_large_1027.mat')

%%
%load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')
%%
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_sbj6_2D.mat')
% traj1 = trj_kxy;
% 
% %%
% load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_sb6.mat')
% traj2 = trj_kxy;
% 
% %%
% a = traj1(:,:,1,1);
% b = traj2(:,:,1,1);
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
%load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleshot_cor_2D_sbj4.mat')
%traj1 = trj_kxy;

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
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{ie} = NUFFTDCF(trj_kxy(:,1,ie), ones(size(DCF.^2)), shift, [120, 120]);
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
Ncoil = 44;%size(data_single, 4); %different

%%
Npar = 24;
Ntime = 175;
img_uncombined = zeros(44,120,120, 6);

%%
data_single = reshape(data_full, [2176, 6, 24, 175, Ncoil]);


%%
atb_2 = zeros(120, 120, 6);
sb = 1;
fprintf('Gridding of k-space data \n')
count = 0;
for ip = 1:Npar
    for it = 1:Ntime
    coil_2 = permute(single(squeeze(coil(:,:,ip,:))), [3,1,2]);
    atb_2 = zeros(120,120,6);
    for ie = 1:6
        DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined(ic,:,:,ie) = par.FT{ie}' * (DCF.*cast(squeeze(data_single(:,ie,ip, it,ic)), par.prec));          
        end
        atb_2(:,:,ie) = sum(img_uncombined(:,:,:,ie).*conj(coil_2), 1);
      
    end
    count = count + 1;
    
    atb_2 = permute(atb_2, [3,1,2]);
    ref_2 = atb_2;
    
%    kspace_2 = kspace_1(:,:,:,:,ii);
%    traj_2 = traj_1;
%    dcf_2 = dcf_1;
%    mask_2 = mask_1;
    idx = count;%+60;
    
    fname = sprintf('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/large_testing/slices_1027_3/slice_%d.mat', count);%+60);
    save(fname, 'atb_2', 'ref_2', 'coil_2', 'idx', 'sb')
    end
end

%%
%img = img_uncombined(:,:,1,1,1);
%figure()
%imshow(abs(img), [])

%%
img_1 = img_uncombined(:,:,:,1,:);
img = permute(img_1, [1,2,5,3,4]);
%%
coil = ZcGetSpiralCoils_2D(img,1);
coil = coil(:,:,:,:,1);

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
for ll = 1:6
    figure()
    echo_num = ll;
for kk = 1:8
    title(['slice', num2str(kk*8-8)])
    subplot(2,4,kk)
    
    imshow(abs(im(:,:,kk*8,echo_num)), [])
    title(['slice', num2str(kk*8)])
end
end

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
