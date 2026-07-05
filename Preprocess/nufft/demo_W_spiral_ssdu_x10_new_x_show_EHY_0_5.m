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

%load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')
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
%%
%Neco = par.Neco;
%Ncoil = par.NcoilFinal;
%Npar = par.Npar;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating N operators \n')
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%%
trj_training = zeros(5, 1088, 6);
trj_validation = zeros(5, 1088, 6);

for m = 1: 5
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(training_all(m,:)~=1) = 0;
    par.FT{m, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [240, 240]);
end
end

for m = 1: 5
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m,:)~=1) = 0;
    par.FT{m+5, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [240, 240]);
end
end

% W_uncombined = zeros([par.imsize*4, Ncoil]);
% M_uncombined = zeros([par.imsize*4, Ncoil]);
% for ip = 1: Npar
%     W_uncombined(:,:,ip) = par.FT{1}' * dcf(:,:,ip);
%     M_uncombined(:,:,ip) = fftshift(fft(ifftshift(W_uncombined(:,:,ip)))) ./ (sqrt((512*2)^2));
% end

M = zeros(10, 6, 240, 240);
for m = 1: 5
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(training_all(m,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(m, ie,:,:) = M1;
end
end

%%
for m = 6: 10
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m-5,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(m, ie,:,:) = M1;
end
end

%%


% W2 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1)))).^0.5;
% M2 = fftshift(fft2(ifftshift(W2))) ./ 240;
% 
% W3 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1))));
% M3 = fftshift(fft2(ifftshift(W3))) ./ 240;

% W1 = par.FT{1}' * (dcR10(:,:,1,1)/max(max(dcR10(:,:,1,1)))).^0;
% M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
%%
% figure()
% imshow(abs(W1), [])
for m = 1: 5
figure()
for ie = 1: 6
    subplot(2,3,ie)
    imshow(squeeze(abs(M(m, ie,:,:))), [])
end
end

%%
figure()
for i = 1:5
    temp = squeeze(M(i,1,:,:));
    subplot(2,3,i)
    imshow(abs(temp), [])

end
% figure()
% imshow(abs(M1), [])
% 
% figure()
% imshow(abs(W2), [])
% 
% figure()
% imshow(abs(M2), [])
% 
% figure()
% imshow(abs(W3), [])
% 
% figure()
% imshow(abs(M3), [])

% %% Gridding reconstruction
% under_mkdata = zeros(512, 600, 4, 6, 44);
% under_mkdata(:, 1:4:end, :, :, :) = mkdata(:, 1:4:end, :, :, :);
% 
% img_uncombined = zeros([par.imsize, Ncoil, Neco, Npar]);
% fprintf('Gridding of k-space data \n')
% for ip = 1:Npar
%     for ie = 1:Neco
%         for ic = 1:Ncoil
%             if par.verbose
%                 fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)
%             end            
%             img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * cast(squeeze(under_mkdata(:,:,ic,ie,ip)), par.prec);          
%         end
%     end
% end
% % img = img_uncombined(:,:,1,1,1);
% % 
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]);
% 
% 
% %%
% under_mkdata_1 = squeeze(under_mkdata(:,:,:,1,:));
% img_uncombined_1 = squeeze(img_uncombined(:,:,:,1,:));



%%
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]

%%
load('/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/LiverData_fromCRettenmeier/dataset_meas_MID00280_FID08424_RAVE.mat')
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')


%%
training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/new_mask/mask_new5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;
%%
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

%%
training_data = zeros(5, 2176, 6, 72, Ncoil);
validation_data  = zeros(5,2176, 6, 72, Ncoil);

%%
for m = 1: 5
    training_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    training_data(m, training_all(m,:)~=1,:,:,:) = 0;
%    validation_data(m, :,:,:,:) = data_single_reshaped(validation_all(m, :)==1, :,:,:);

    validation_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
    validation_data(m, validation_all(m,:)~=1,:,:,:) = 0;
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
%            img_uncombined_validation(m, :,:,ic,ie,ip) = par.FT{m+5, ie}' * (DCF.*cast(squeeze(validation_data(m, :,ie,ip,ic)), par.prec));
        end
    end
end
end
%%
echo1_train = squeeze(img_uncombined_train(:,:,:,:,1,:));
echo1_validation = squeeze(img_uncombined_validation(:,:,:,:,1,:));

%%
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb7_coil.mat')

%%
coil_original = coil;
coil = permute(coil, [1,2,4,3]);
sense1 = zeros(5,120,120,72);
for i = 1: 5
    sense1(i,:,:,:) = sum(squeeze(echo1_train(i,:,:,:,:)).*conj(coil), 3);
end

%%
cg_rec = zeros(5,120,120,72);
for i = 1:5
    cg_rec(i,:,:,:) = ZcCGSENSE3D_SingleEcho(squeeze(sense1(i,:,:,:)), coil_original,0,squeeze(M(i,1,:,:))/240,10);
end

%%
figure()
imshow(abs(squeeze(cg_rec(1,:,:,34))), [])

figure()
imshow(abs(squeeze(sense1(1,:,:,34))), [])

diff = squeeze(sense1(1,:,:,34)) - cg_rec(1,:,:,34);

a = max(max(abs(sense1(1,:,:,34))));
b = max(max(abs(cg_rec(1,:,:,34))));

sense1_28 = permute(sense1(:,:,:,28), [2,3,1]);
sense1_34 = permute(sense1(:,:,:,34), [2,3,1]);

cg_28 = permute(cg_rec(:,:,:,28), [2,3,1]);
cg_34 = permute(cg_rec(:,:,:,34), [2,3,1]);

sense1_28 = reshape(sense1_28, [120,600]);
sense1_34 = reshape(sense1_34, [120,600]);
cg_28 = reshape(cg_28, [120,600]);
cg_34 = reshape(cg_34, [120,600]);

%%
figure()
imshow(abs(sense1_28), [])

figure()
imshow(abs(cg_28), [])

figure()
imshow(abs(sense1_34), [])

figure()
imshow(abs(cg_34), [])

%%
M_e1 = squeeze(M(:,1,:,:));
M_e1 = permute(M_e1, [2,3,1]);

M_e1 = reshape(M_e1,[240, 240*5]);

figure()
imshow(abs(M_e1), [])

%%
a = [sense1_28; cg_28; sense1_34; cg_34];

figure()
imshow(abs(a), [])

%%
b = [cg_28, cg_34];
figure()
imshow(abs(b), [])

%%
coil_34 = squeeze(coil(:,:,34,:));

img_uncombined_34 = echo1_train(:,:,:,:,34);

sense1 = zeros(120,120,5);

for i = 1:5
    sense1(:,:,i) = sum(squeeze(img_uncombined_34(i,:,:,:)).*conj(coil_34), 3);
end

%figure()
reshaped_sense1_train_34 = reshape(sense1, [120, 600]);

%imshow(abs(reshaped_sense1), [])
%%
img_uncombined_34_test = echo1_validation(:,:,:,:,34);

sense1 = zeros(120,120,5);

for i = 1:5
    sense1(:,:,i) = sum(squeeze(img_uncombined_34_test(i,:,:,:)).*conj(coil_34), 3);
end

%figure()
reshaped_sense1_test_34 = reshape(sense1, [120, 600]);

%%
coil_28 = squeeze(coil(:,:,28,:));

img_uncombined_28 = echo1_train(:,:,:,:,28);

sense1 = zeros(120,120,5);

for i = 1:5
    sense1(:,:,i) = sum(squeeze(img_uncombined_28(i,:,:,:)).*conj(coil_28), 3);
end

%figure()
reshaped_sense1_train_28 = reshape(sense1, [120, 600]);

%imshow(abs(reshaped_sense1), [])
%%
img_uncombined_28_test = echo1_validation(:,:,:,:,28);

sense1 = zeros(120,120,5);

for i = 1:5
    sense1(:,:,i) = sum(squeeze(img_uncombined_28_test(i,:,:,:)).*conj(coil_28), 3);
end

%figure()
reshaped_sense1_test_28 = reshape(sense1, [120, 600]);

%%
figure()
im = [reshaped_sense1_train_34; reshaped_sense1_test_34;reshaped_sense1_train_28; reshaped_sense1_test_28];
imshow(abs(im), [])
%%
figure()
imshow(abs(reshaped_sense1_train_34) ,[])

%%
load /home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/training_slices_new/slice_1
%%
%coil_2 = permute(coil_34, [3,2,1]);

%coil_2 = permute(coil_28, [3,2,1]);

function  u = ZcCGSENSE3D_SingleEcho(FHy, coil, alpha,RapidOperator,numIter)

%y = FHy/max(abs(FHy(:)));
y = FHy;
[nx,ny,nz,nc] = size(coil);
M  = @(x) applyM(RapidOperator,coil,x) + alpha*x;
x = 0*y(:);
r = y(:);
p = r;
rr = r'*r;

list = [FHy(:,:,40)];
previous = FHy(:,:,40);
difflist = zeros([120,120],'single');
%L2list = [];

for it = 1:numIter
    Ap = M(p);
    a = rr/(p'*Ap);
    x = x + a*p;
    %L2list = [L2list, (r'*r) / (a'*a);];
    %aa=reshape(a*p,nx,ny,nz);
    %figure;imshow(abs(squeeze(aa(:,64,:))),[])
    rnew = r - a*Ap;
    b = (rnew'*rnew)/rr;
    r=rnew;
    rr = r'*r;
    p = r + b*p;
    %disp([num2str(max(abs(col(a)))),';',num2str(max(abs(col(r)))),';',num2str(max(abs(col(b)))),';',num2str(max(abs(col(p)))),';']);

    u_it = reshape(x,nx,ny,nz);
    diff = abs(u_it(:,:,40) - previous)*5;
    list = [list,u_it(:,:,40)];
    difflist = [difflist, diff];
    previous = u_it(:,:,40);
end
%figure;imshow(abs([list;difflist]),[])
u  = reshape(x,nx,ny,nz);
%figure;plot(L2list)
end

%% Derivative evaluation
function y = applyM(M,coil,x)
[nx,ny,slice,nc] = size(coil);
dx = reshape(x,[nx,ny,slice,1]);
ZeroPadded_img = zeros([nx*2,ny*2,slice,nc],'single');
coilimg = repmat(dx,[1,1,1,nc]).*coil;
ZeroPadded_img(1+nx/2:nx/2+nx,1+ny/2:ny/2+ny,:,:)=coilimg;
%size(M)
%size(ZeroPadded_img)
itermediate = repmat(M,[1,1,slice,nc]).*  fft(fft(ZeroPadded_img ,[],1),[],2)/sqrt(2*nx*2*ny);
rapid  = ifft(ifft( itermediate ,[],1),[],2) * sqrt(2*nx*2*ny);
y = (sum(rapid(1+nx/2:nx/2+nx,1+ny/2:ny/2+ny,:,:).*conj(coil),4));
%figure;imshow(abs(y),[]);title('cg intermediate')
y = y(:) ;
end
