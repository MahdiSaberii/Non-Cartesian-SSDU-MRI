%% Singleshot data
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')

%% DCF
load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF.mat')

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating NUFFT operators \n')

for ie = 1:Neco
    FT{ie} = NUFFTDCF(trj_kxy(:,1,ie), ones(size(trj_kxy(:,1,ie))), shift, [120, 120]);
end

Ncoil = size(data_single, 3);
Npar = size(data_single, 4);
img_uncombined = zeros([[120,120], Ncoil, Neco, Npar]);
data_single_reshaped = reshape(data_single, [2176, 6, 72, Ncoil]);

%% y
y = squeeze(data_single_reshaped(:,1,40,:));

%% my_cg_out
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/training_slices_simple/cg_10x_example3_M_1.mat' %my cg-out
my_recon = squeeze(recon(1,1,:,:));

%% load M
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/M_10x_1.mat' %M_1 represents F_2N*delta

%% My cg input
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/testing_slices_10x_example3/slice_1.mat' %my cg input

%% z (zero-filled) w_W
z_w_W = squeeze(atb_2(1,:,:));

%% Chi's cg
slice_coil = zeros(120,120,1,44);
slice_coil(:,:,1,:) = permute(coil_2,[2,3,1]);
chi_recon = ZcCGSENSE3D_SingleEcho(permute(atb_2(1,:,:), [2,3,1]), slice_coil, 0,squeeze(M(1,:,:)),10); % chi's cg-out

%% plot sense1 and my_cg together
atb = squeeze(abs(atb_2(1,:,:)));
figure()
imshow(abs([atb, my_recon]), [])

%% Ex
idx = 40;

coil_part = permute(coil_2, [2,3,1]);
coil_img = coil_part.*my_recon;
Ex_my_recon = zeros(2176, 44);
for i = 1:Ncoil
    Ex_my_recon(:,i) = FT{1}*coil_img(:,:,i);
end

%%
label = squeeze(data_single_reshaped(:,1,idx,:));

%% F
FT_1 = FT{1};

%% EHEx_o_w
uncombined_EHE_o_W = zeros(120, 120, 44);

for i = 1: Ncoil
    uncombined_EHE_o_W(:,:,i) = FT{1}'*Ex_my_recon(:,i);
end
EHE_o_W = sum(uncombined_EHE_o_W.*conj(coil_part), 3);

%% EHEx_w_W
uncombined_EHE_w_W = zeros(120, 120, 44);
DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
for i = 1: Ncoil
    uncombined_EHE_w_W(:,:,i) = FT{1}'*(DCF.*Ex_my_recon(:,i));
end
EHE_w_W = sum(uncombined_EHE_w_W.*conj(coil_part), 3);


%% EHEx_fast
M_e1 = squeeze(M(1,:,:));

coil_part = permute(coil_2, [2,3,1]);
coil_img = coil_part.*my_recon;

padded = zeros(240,240,44);
padded(61:180, 61:180, :) = coil_img;

fft_domain_multiplied_by_M = zeros(240,240,44);
ifft_domain = zeros(240,240,44);
for i = 1: 44
    fft_domain_multiplied_by_M(:,:,i) = fftshift(fft2(padded(:,:,i))) ./ 240 .*M_e1;
    ifft_domain(:,:,i) = ifftshift(ifft2(fft_domain_multiplied_by_M(:,:,i))) .* 240;
end
truncated = ifft_domain(61:180, 61:180, :);
EHE_fast = sum(uncombined_EHE_w_W.*conj(truncated), 3);


%%
idx = 40;
coil_part = permute(coil_2, [2,3,1]);

coil_img = coil_part.*chi_recon;
ksp_vec_chi_recon = zeros(2176, 44);
for i = 1:Ncoil
    ksp_vec_chi_recon(:,i) = FT{1}*coil_img(:,:,i);
end

%%
DCF = dc10x(:,1)./ max(max(dc10x(:,1)));
fourier = zeros(120,120,Ncoil);
for i = 1:44
    fourier(:,:,i) = FT{1}' * (DCF.*label(:,i));
end

sense1 = sum(fourier.*conj(coil_part), 3);
max_sense1 = max(max(abs(sense1)));

%% z_o_W
label = squeeze(data_single_reshaped(:,1,40,:));

DCF = dc10x(:,1)./ max(max(dc10x(:,1)));
fourier = zeros(120,120,Ncoil);
for i = 1:44
    fourier(:,:,i) = FT{1}' * label(:,i);
end

z_o_W = sum(fourier.*conj(coil_part), 3);
%max_sense1 = max(max(abs(sense1)));


%%
max_sense1 = max(max(abs(atb_2(1,:,:))));

%%
my_max = max(max(abs(my_recon)));
chi_max = max(max(abs(chi_recon)));

%%
ksp_vec = ksp_vec_chi_recon;


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
