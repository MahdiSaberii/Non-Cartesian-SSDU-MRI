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
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')

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

load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb7_coil.mat')
Ncoil = size(data_single, 3); %different
Npar = 72;
img_uncombined = zeros([[120,120], Ncoil, 6, Npar]);
data_single_reshaped = reshape(data_single, [2176, 6, 72, Ncoil]);

%%
%load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/training_slices_simple/cg_10x_example2_M_1.mat'

load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/M_10x_1.mat'
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/testing_slices_10x_example2/slice_1.mat'
coil_temp = zeros(120,120,1,44);
coil_temp(:,:,1,:) = permute(coil_2,[2,3,1]);
recon = ZcCGSENSE3D_SingleEcho(permute(atb_2(1,:,:), [2,3,1]), coil_temp, 0,squeeze(M(1,:,:)),10);


%%
idx = 40;
coil_part = squeeze(coil(:,:,idx,:));

coil_img1 = recon.*coil_part;
ksp_vec = zeros(2176, 44);
for i = 1:44
    ksp_vec(:,i) = par.FT{1}*coil_img1(:,:,i);
end

label = squeeze(data_single_reshaped(:,1,idx,:));

%%
DCF = dc10x(:,1)./ max(max(dc10x(:,1)))*2;
fourier = zeros(120,120,44);
for i = 1:44
    fourier(:,:,i) = par.FT{1}' * (DCF.*label(:,i));
end

sense1 = sum(fourier.*coil_part, 3);
a = max(max(abs(sense1)));

%%
recon1 = squeeze(recon(1,1,:,:));
b = max(max(abs(recon1)));

c = max(max(abs(par.FT{1}' * (DCF.*(par.FT{1} * recon1)))));
d = max(max(abs(par.FT{1}' * (par.FT{1} * recon1))));


%%
diff1 = ksp_vec-label;

%%
e = max(max(abs(ksp_vec)));
f = max(max(abs(label)));

%%
DCF = dc10x(:,1)./ max(max(dc10x(:,1)));

diff2 = diff1.*sqrt(DCF);

normalization1 = norm(diff1(:))/norm(label(:));
normalization2 = norm(diff1(:), 1)/norm(label(:), 1);

label1 = label.*sqrt(DCF);

normalization3 = norm(diff2(:))/norm(label1(:));
normalization4 = norm(diff2(:), 1)/norm(label1(:), 1);

%%
factor = mean(abs(label(:)))./mean(abs(ksp_vec(:)));

%%
ksp_vec_new = ksp_vec.*factor;
diff1 = ksp_vec_new-label;
diff2 = diff1.*sqrt(DCF);

normalization5 = norm(diff1(:))/norm(label(:));
normalization6 = norm(diff1(:), 1)/norm(label(:), 1);

label1 = label.*sqrt(DCF);

normalization7 = norm(diff2(:))/norm(label1(:));
normalization8 = norm(diff2(:), 1)/norm(label1(:), 1);

%%
a = [ksp_vec, label];

figure()
imshow(abs(a.'), [])

%%
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:Neco
        for ic = 1:Ncoil
            fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)           
            img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * cast(squeeze(data_single_reshaped(:,ie,ip,ic)), par.prec);          
        end
    end
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

%%
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]
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

function  u = ZcCGSENSE3D_SingleEcho(FHy, coil, alpha,RapidOperator,numIter)

%y = FHy/max(abs(FHy(:)));
y = FHy;
[nx,ny,nz,nc] = size(coil);
M  = @(x) applyM(RapidOperator,coil,x) + alpha*x;
x = 0*y(:);
r = y(:);
p = r;
rr = r'*r;
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
