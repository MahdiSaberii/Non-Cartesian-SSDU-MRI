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

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
% a = dcR10(:,1,1,1);
% figure()
% plot(a)


%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/full_sb7.mat')
%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/1_iter_cg.mat')

%%
data_full_reshaped = reshape(data_full, [2176, 6, 10, 72, 44]);

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb7_coil.mat')
%%
shift   = [0,0];
trj = trj_kxy(:,:,1,1);
trj_part = trj_kxy(:,1:5:end,1,1);
par.FT{1} = NUFFTDCF(trj, ones(size(trj)), shift, [120, 120]);

par.FT{2} = NUFFTDCF(trj_part, ones(size(trj_part)), shift, [120, 120]);

%%
input_full = squeeze(data_full_reshaped(:,1,:,34,1));
input_part = squeeze(data_full_reshaped(:,1,1:5:end,34,1));
img_full_o = par.FT{1}' * input_full;
img_full_w = par.FT{1}' * (dcFull(:,:,1).*input_full);

img_part_o = par.FT{2}' * input_part;
img_part_w1 = par.FT{2}' * (dcR5(:,:,1, 1).*input_part);
img_part_w2 = par.FT{2}' * (dcFull(:,1:5:end,1).*input_part);

figure()
subplot(2,3,1)
imshow(abs(img_full_o),[])

subplot(2,3,2)
imshow(abs(img_full_w),[])

subplot(2,3,3)
imshow(abs(img_part_o),[])

subplot(2,3,4)
imshow(abs(img_part_w1),[])

subplot(2,3,5)
imshow(abs(img_part_w2),[])

normalization1 = norm(img_full_o).^2;
normalization2 = norm(img_full_w).^2;
normalization3 = norm(img_part_o).^2;
normalization4 = norm(img_part_w1).^2;
normalization5 = norm(img_part_w2).^2;

max1 = max(max(abs(img_full_o)));
max2 = max(max(abs(img_full_w)));
max3 = max(max(abs(img_part_o)));
max4 = max(max(abs(img_part_w1)));
max5 = max(max(abs(img_part_w2)));

%%
coil_part = squeeze(coil(:,:,34,:));

coil_img1 = squeeze(cg_sense1(1,1,:,:)).*coil_part;
coil_img3 = squeeze(cg_sense3(1,1,:,:)).*coil_part;

ksp_vec1 = zeros(2176, 2, 44);

ksp_vec2 = zeros(2176, 10, 44);

ksp_vec3 = zeros(2176, 10, 44);

for i = 1:44
    ksp_vec1(:,:,i) = par.FT{2}*coil_img1(:,:,i);
    ksp_vec3(:,:,i) = par.FT{1}*coil_img3(:,:,i);
end

DCF_1 = dcR5(:,:,1,1)./max(max(dcR5(:,:,1,1)));
DCF_2 = dcFull(:,:,1)./max(max(dcFull(:,:,1)));

label_ksp1 = squeeze(data_full_reshaped(:,1,1:5:end,34,:));
label_ksp2 = squeeze(data_full_reshaped(:,1,:,34,:));

diff1 = ksp_vec1-label_ksp1;
diff2 = ksp_vec2-label_ksp2;

diff1_temp = diff1.*sqrt(DCF_1);
diff2_temp = diff2.*sqrt(DCF_2);

normalization1 = norm(diff1(:)).^2;
normalization1_1 = norm(diff1_temp(:)).^2;

normalization2 = norm(diff2(:)).^2;
normalization2_2 = norm(diff2_temp(:)).^2;

%%
dc_f = dcFull(:,:,1)./max(max(abs(dcFull(:,:,1))));
dc_fp = dcFull(:,1:5:end,1)./max(max(abs(dcFull(:,1:5:end,1))));
dc_p = dcR5(:,:,1, 1)./ max(max(abs(dcR5(:,:,1, 1))));

%%

input_full = squeeze(data_full_reshaped(:,1,:,34,1));
input_part = squeeze(data_full_reshaped(:,1,1:5:end,34,1));
img_full_o = par.FT{1}' * input_full;
img_full_w = par.FT{1}' * (dc_f.*input_full);

img_part_o = par.FT{2}' * input_part;
img_part_w1 = par.FT{2}' * (dc_p.*input_part);
img_part_w2 = par.FT{2}' * (dc_fp.*input_part);

figure()
subplot(2,3,1)
imshow(abs(img_full_o),[])

subplot(2,3,2)
imshow(abs(img_full_w),[])

subplot(2,3,3)
imshow(abs(img_part_o),[])

subplot(2,3,4)
imshow(abs(img_part_w1),[])

subplot(2,3,5)
imshow(abs(img_part_w2),[])

normalization1 = norm(img_full_o).^2;
normalization2 = norm(img_full_w).^2;
normalization3 = norm(img_part_o).^2;
normalization4 = norm(img_part_w1).^2;
normalization5 = norm(img_part_w2).^2;

max1 = max(max(abs(img_full_o)));
max2 = max(max(abs(img_full_w)));
max3 = max(max(abs(img_part_o)));
max4 = max(max(abs(img_part_w1)));
max5 = max(max(abs(img_part_w2)));

%%


% %%
% figure()
% a = imag(trj_kxy(:, 1, 1));
% plot(a)
% shift   = [0,0];
% vec = zeros(2176,1);
% vec(1,1) = 1;
% 
% traj = trj_kxy(:,1,:);
% 
% for ie = 1:6
%     par.FT{ie} = NUFFTDCF(traj(:,:,ie), ones(size((dc10x(:,ie)./max(max(dc10x(:,ie)))).^2)), shift, [240, 240]);
% end
% 
% result = par.FT{1}' * vec;
% 
% delta_mat1 = zeros(240, 240);
% %delta_mat1(121,121) = 1;
% delta_mat1(121,121) = 240;
% 
% example1 = par.FT{1} * delta_mat1;
% 
% normalization1 = norm(example1).^2;
% 
% DCF1 = dc10x(:,1)./ max(max(dc10x(:,1)));
% DCF2 = dc10x(:,1);
% 
% example2 = par.FT{1}' * (DCF1.*example1);
% example3 = par.FT{1}' * (DCF2.*example1);
% 
% normalization2 = norm(example2).^2;
% normalization3 = norm(example3).^2;
% 
% example4 = par.FT{1}' * example1;
% normalization4 = norm(example4).^2;
% 
% %%
% max4 = max(max(abs(example4)));
% max2 = max(max(abs(example2)));
% 
% % %%
% % figure()
% % imshow(squeeze(abs(result)), [])
% % 
% % figure()
% % imshow(squeeze(abs(result)), [0.0041 0.0042])
% % 
% % %%
% % %Neco = par.Neco;
% % %Ncoil = par.NcoilFinal;
% % %Npar = par.Npar;
% % 
% % delta_mat1 = zeros(240, 240);
% % delta_mat1(121,121) = 1;
% % 
% % delta_mat2 = zeros(240, 240);
% % delta_mat2(120,120) = 1;
% % 
% % delta_mat3 = zeros(240, 240);
% % delta_mat3(1, 1) = 1;
% % 
% % two_shot_traj = squeeze(trj_kxy(:,1:5:end,1,:));
% % 
% % 
% % %% Calculate NUFFT operator
% % shift   = [0,0];
% % Neco = 6;
% % fprintf('Calculating NUFFT operators \n')
% % 
% % for ie = 1:Neco
% %     par.FT{ie} = NUFFTDCF(two_shot_traj(:,:,ie), ones(size((dcR5(:,:,ie,1)./max(max(dcR5(:,:,ie,1)))).^2)), shift, [240, 240]);
% % end
% % 
% % example1 = par.FT{ie} * delta_mat1;
% % std_delta1 = std(example1);
% % 
% % example2 = par.FT{ie} * delta_mat2;
% % std_delta2 = std(example2);
% % 
% % example3 = par.FT{ie} * delta_mat3;
% % std_delta3 = std(example3);
% % 
% % mean_delta = mean(mean(example1));
% % normalized_delta = zeros(240, 240);
% % normalized_delta(121,121) = 1000/mean_delta;
% % 
% % example4 = par.FT{ie} * normalized_delta;
% % std_delta4 = std(example4);
% % 
% % % W_uncombined = zeros([par.imsize*4, Ncoil]);
% % % M_uncombined = zeros([par.imsize*4, Ncoil]);
% % % for ip = 1: Npar
% % %     W_uncombined(:,:,ip) = par.FT{1}' * dcf(:,:,ip);
% % %     M_uncombined(:,:,ip) = fftshift(fft(ifftshift(W_uncombined(:,:,ip)))) ./ (sqrt((512*2)^2));
% % % end
% % 
% % M = zeros(6, 240, 240);
% % for ie = 1:Neco
% %     DCF = dcR5(:,:,ie,1)./ max(max(dcR5(:,:,ie,1)));
% %     
% %     intermediate = DCF.* (par.FT{ie} * (delta_mat1*240));
% %     
% %     W1 = par.FT{ie}' * intermediate;
% %     M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
% %     M(ie,:,:) = M1;
% % end
% % 
% % % W2 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1)))).^0.5;
% % % M2 = fftshift(fft2(ifftshift(W2))) ./ 240;
% % % 
% % % W3 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1))));
% % % M3 = fftshift(fft2(ifftshift(W3))) ./ 240;
% % 
% % % W1 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1)))).^0;
% % % M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
% % %%
% % % figure()
% % % imshow(abs(W1), [])
% % figure()
% % for ie = 1: 6
% %     subplot(2,3,ie)
% %     imshow(squeeze(abs(M(ie,:,:))), [0 1e-1])
% % end
% % 
% % figure()
% % for ie = 1: 6
% %     subplot(2,3,ie)
% %     imshow(squeeze(angle(M(ie,:,:))), [-pi, pi])
% % end
% % % figure()
% % % imshow(abs(M1), [])
% % % 
% % % figure()
% % % imshow(abs(W2), [])
% % % 
% % % figure()
% % % imshow(abs(M2), [])
% % % 
% % % figure()
% % % imshow(abs(W3), [])
% % % 
% % figure()
% % imshow(abs(M3), [])
% 
% % %% Gridding reconstruction
% % under_mkdata = zeros(512, 600, 4, 6, 44);
% % under_mkdata(:, 1:4:end, :, :, :) = mkdata(:, 1:4:end, :, :, :);
% % 
% % img_uncombined = zeros([par.imsize, Ncoil, Neco, Npar]);
% % fprintf('Gridding of k-space data \n')
% % for ip = 1:Npar
% %     for ie = 1:Neco
% %         for ic = 1:Ncoil
% %             if par.verbose
% %                 fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)
% %             end            
% %             img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * cast(squeeze(under_mkdata(:,:,ic,ie,ip)), par.prec);          
% %         end
% %     end
% % end
% % % img = img_uncombined(:,:,1,1,1);
% % % 
% % img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% % im = permute(squeeze(img_gridding_sos),[1 2 4 3]);
% % 
% % 
% % %%
% % under_mkdata_1 = squeeze(under_mkdata(:,:,:,1,:));
% % img_uncombined_1 = squeeze(img_uncombined(:,:,:,1,:));
% 
% 
% 
% %%
% % img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% % im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]
% aa = 1:(2176*2);
% figure()
% plot(aa, transpose(real(example3(:))))
% 
% figure()
% plot(aa, transpose(imag(example3(:))))
