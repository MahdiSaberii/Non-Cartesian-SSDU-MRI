load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/sb6_coil_new_data_modified_320.mat')

%%
part = img_uncombined(:,:,:,1,34);

img = zeros(480, 1560);
for i = 1:4
    for j = 1:13
        img(1+(i-1)*120: i*120, 1+(j-1)*120:j*120) = part(:,:,(i-1)*10+j);
    end
end

figure()
imshow(abs(img), [])

%%
ksp = zeros(size(part));
for i = 1:52
    ksp(:,:,i) = fftshift(fft2(ifftshift(part(:,:,i))));
end

kspace = zeros(480, 1560);
for i = 1:4
    for j = 1:13
        kspace(1+(i-1)*120: i*120, 1+(j-1)*120:j*120) = ksp(:,:,(i-1)*10+j);
    end
end

figure()
imshow(log(abs(kspace)), [])

%%
figure()
imshow(angle(kspace), [])

%%
load('/home/daedalus1-raid1/chi/Hawaii/PGDL/database/Testing/Coils_2D/sb6/sb_6_slice_34.mat')

img = zeros(480, 1560);
for i = 1:4
    for j = 1:13
        img(1+(i-1)*120: i*120, 1+(j-1)*120:j*120) = coil(:,:,(i-1)*10+j);
    end
end

figure()
imshow(abs(img), [])
%%
part = squeeze(coil(:,:,1,:));

img = zeros(320*4, 320*13);
for i = 1:4
    for j = 1:13
        img(1+(i-1)*320: i*320, 1+(j-1)*320:j*320) = part(:,:,(i-1)*10+j);
    end
end

figure()
imshow(abs(img), [])