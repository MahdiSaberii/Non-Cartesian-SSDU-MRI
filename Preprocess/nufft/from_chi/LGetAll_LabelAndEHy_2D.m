addpath('/home/daedalus1-raid1/chi/kooshball/gpuNUFFT-master/CUDA/bin');
addpath(genpath('/home/daedalus1-raid1/chi/kooshball/gpuNUFFT-master/matlab/demo/utils'));
addpath(genpath('/home/daedalus1-raid1/chi/kooshball/gpuNUFFT-master/gpuNUFFT'));
addpath(genpath('ESPIRiT_original/'));

osf = 2; % oversampling: 1.5 1.25
wg = 5; % kernel width: 5 7 
sw = 8; % parallel sectors' width: 12 16

load('../DCF.mat')
for ii = 1:size(dcFull,3)
    dcFull(:,:,ii) = dcFull(:,:,ii)/max(max(dcFull(:,:,ii)));
end
dcFull = reshape(dcFull,[2176*10,6]);
load('../traj_full.mat')
imwidth = 120;


for sbj = 1:6
    % Firstly determine the sbj to be used in trainig or testing
    if sbj<6
        Category = 'Training'
    else
        Category = 'Testing' 
    end

    disp(['Loading Subject ',num2str(sbj)]);
    load(['../RawData2D/full_sbj',num2str(sbj),'_2D.mat']);
    data_full = reshape(data_full, [2176, 6, size(data_full, 2) size(data_full, 3), size(data_full, 4)]) ;
    disp(['---------> Loading is done. Firstly get gridding... <----------'])
    [Nspoke, Necho, Nleaves, Nslice, Ncoil] = size(data_full);
    img = zeros([imwidth,imwidth,Nslice,Ncoil,Necho],'single');
    ksp = reshape(permute(data_full,[1,3,2,4,5]),[Nspoke*Nleaves,Necho,Nslice,Ncoil]);
    
    for IdxEcho = 1:Necho
        trajx = kx(:,:,IdxEcho);
        trajy = ky(:,:,IdxEcho);
        FT = gpuNUFFT([trajx(:),trajy(:)]',ones(size(dcFull(:,IdxEcho))),osf,wg,sw,[imwidth,imwidth],[]);
        for idxSlice = 1:Nslice
            for idxCoil = 1:Ncoil

                img(:,:,idxSlice,idxCoil,IdxEcho) = FT'*(dcFull(:,IdxEcho).*ksp(:,IdxEcho,idxSlice,idxCoil));
                disp(['Echo #',num2str(IdxEcho),', Coil #',num2str(idxCoil),', Slice #', num2str(idxSlice)])
            end
        end
    end
    
    disp(['---------> Get coil maps ... <----------'])
    coil = ZcGetSpiralCoils_2D(img,1);
    coil = coil(:,:,:,:,1);
    for IdxCoil = 1:size(coil,4)
        fileName = ['../',Category,'/Coils_2D/sb',num2str(sbj),'/sb_',num2str(sbj),'_coil_',num2str(IdxCoil),'.mat']
        coil_singleChannel = single(coil(:,:,:,IdxCoil));
        save(fileName,'coil_singleChannel');
    end
    
    disp(['---------> Get label images <----------'])
    label = sum(img.*repmat(conj(coil),[1,1,1,1,Necho]),4);
    
    % Use CGSENSE with M to adjust scaling using M. 
    load('../FastOperator/Full.mat')
    for IdxEcho = 1:Necho
        label(:,:,:,:,IdxEcho) = ZcCGSENSE3D_SingleEcho(label(:,:,:,:,IdxEcho), coil, 0,M_Full(:,:,IdxEcho),1);
    end
    % Save the data
    label = cat(4,real(label),imag(label)); % X, Y, Slice, R/I, Echo
    %label = permute(label,[5,4,3,2,1]); % Echo, R/I, Slice, Y, X
    for IdxSlice = 1:Nslice
        label_SingleSlice = single(squeeze(label(:,:,IdxSlice,:,:)));
        fileName = ['../',Category,'/Label_2D/sb',num2str(sbj),'/sb_',num2str(sbj),'_slice_',num2str(IdxSlice),'.mat']
        save(fileName,'label_SingleSlice');
        for IdxEcho = 1:Necho
            slice_img = squeeze(label_SingleSlice(:,:,1,IdxEcho) + label_SingleSlice(:,:,2,IdxEcho)*1j);
            slice_ksp = log(abs(fftshift(fft2(slice_img))));
            slice_ksp = slice_ksp - min(slice_ksp(:));
            slice_ksp = slice_ksp/max(slice_ksp(:));
            slice_img = abs(slice_img);
            slice_img = slice_img/max(slice_img(:));
            PNG2Bsaved = [slice_img,slice_ksp];
            PNG2Bsaved = PNG2Bsaved/max(PNG2Bsaved(:));
            PNGName = ['../',Category,'/Label_2D/sb',num2str(sbj),'/PNG/sb_',num2str(sbj),'_slice_',num2str(IdxSlice),'_Echo_',num2str(IdxEcho),'.png'];
            imwrite(PNG2Bsaved,PNGName);
        end
    end
    
    disp('---------> Get R = 5 data <----------');
    load('../FastOperator/R5.mat');
    Nleaves = 2;
    for Traj_Set = 1:5
        img = zeros([imwidth,imwidth,Nslice,Ncoil,Necho],'single');
        ksp = reshape(permute(data_full(:,:,Traj_Set:5:end,:,:),[1,3,2,4,5]),[Nspoke*Nleaves,Necho,Nslice,Ncoil]);

        for IdxEcho = 1:Necho
            dcf = dcR5(:,:,IdxEcho,Traj_Set);
            dcf = dcf/max(dcf(:));
            dcf = dcf(:);        
            trajx = kx(:,Traj_Set:5:end,IdxEcho);
            trajy = ky(:,Traj_Set:5:end,IdxEcho);
            FT = gpuNUFFT([trajx(:),trajy(:)]',ones(size(dcf)),osf,wg,sw,[imwidth,imwidth],[]);
            for idxSlice = 1:Nslice
                for idxCoil = 1:Ncoil

                    img(:,:,idxSlice,idxCoil,IdxEcho) = FT'*(dcf.*ksp(:,IdxEcho,idxSlice,idxCoil));
                    disp(['Echo #',num2str(IdxEcho),', Coil #',num2str(idxCoil),', Slice #', num2str(idxSlice)])
                end
            end
        end
        Ehy = sum(img.*repmat(conj(coil),[1,1,1,1,Necho]),4);
        cgsense = zeros(size(Ehy),'single');
        % Run CGSENSE, 5 Iters, fix the scaling & a little bit dealiasing
        M = M_R5(:,:,:,Traj_Set);
        for IdxEcho = 1:Necho
            cgsense(:,:,:,:,IdxEcho) = ZcCGSENSE3D_SingleEcho(Ehy(:,:,:,:,IdxEcho), coil, 0.01, M(:,:,IdxEcho),5);
        end
        Ehy = cat(4,real(Ehy),imag(Ehy)); % X, Y, Slice, R/I, Echo
        %Ehy = permute(Ehy,[5,4,3,2,1]); % Echo, R/I, Slice, Y, X
        cgsense = cat(4,real(cgsense),imag(cgsense)); % X, Y, Slice, R/I, Echo
        %cgsense = permute(cgsense,[5,4,3,2,1]); % Echo, R/I, Slice, Y, X
        
        for IdxSlice = 1:Nslice
            EHy_SingleSlice = single(squeeze(Ehy(:,:,IdxSlice,:,:)));
            cgsense_SingleSlice = single(squeeze(cgsense(:,:,IdxSlice,:,:)));
            fileName = ['../',Category,'/EHy_R5_2D/sb',num2str(sbj),'/traj_',num2str(Traj_Set),'_sb_',num2str(sbj),'_slice_',num2str(IdxSlice),'.mat']
            save(fileName,'EHy_SingleSlice', 'cgsense_SingleSlice');
            for IdxEcho = 1:Necho
                slice_img = squeeze(EHy_SingleSlice(:,:,1,IdxEcho) + EHy_SingleSlice(:,:,2,IdxEcho)*1j);
                slice_ksp = log(abs(fftshift(fft2(slice_img))));
                slice_ksp = slice_ksp - min(slice_ksp(:));
                slice_ksp = slice_ksp/max(slice_ksp(:));
                slice_img = abs(slice_img);
                slice_img = slice_img/max(slice_img(:));
                PNG2Bsaved = [slice_img,slice_ksp];
                PNG2Bsaved = PNG2Bsaved/max(PNG2Bsaved(:));
                PNGName = ['../',Category,'/EHy_R5_2D/sb',num2str(sbj),'/PNG/traj_',num2str(Traj_Set),'_sb_',num2str(sbj),'_slice_',num2str(IdxSlice),'_Echo_',num2str(IdxEcho),'.png'];
                imwrite(PNG2Bsaved,PNGName);
            end
        end
    end
    
    disp('---------> Get R = 10 data <----------');
    load('../FastOperator/R10.mat');
    Nleaves = 1;
    for Traj_Set = 1:10
        img = zeros([imwidth,imwidth,Nslice,Ncoil,Necho],'single');
        ksp = reshape(permute(data_full(:,:,Traj_Set:10:end,:,:),[1,3,2,4,5]),[Nspoke*Nleaves,Necho,Nslice,Ncoil]);

        for IdxEcho = 1:Necho
            dcf = dcR10(:,:,IdxEcho,Traj_Set);
            dcf = dcf/max(dcf(:));
            dcf = dcf(:);
            trajx = kx(:,Traj_Set:10:end,IdxEcho);
            trajy = ky(:,Traj_Set:10:end,IdxEcho);
            FT = gpuNUFFT([trajx(:),trajy(:)]',ones(size(dcf)),osf,wg,sw,[imwidth,imwidth],[]);
            for idxSlice = 1:Nslice
                for idxCoil = 1:Ncoil

                    img(:,:,idxSlice,idxCoil,IdxEcho) = FT'*(dcf.*ksp(:,IdxEcho,idxSlice,idxCoil));
                    disp(['Echo #',num2str(IdxEcho),', Coil #',num2str(idxCoil),', Slice #', num2str(idxSlice)])
                end
            end
        end
        Ehy = sum(img.*repmat(conj(coil),[1,1,1,1,Necho]),4);
        cgsense = zeros(size(Ehy),'single');
        
        % Run CGSENSE, 5 Iters, fix the scaling & a little bit dealiasing
        M = M_R10(:,:,:,Traj_Set);
        for IdxEcho = 1:Necho
            cgsense(:,:,:,:,IdxEcho) = ZcCGSENSE3D_SingleEcho(Ehy(:,:,:,:,IdxEcho), coil, 0.01, M(:,:,IdxEcho),5);
        end
        Ehy = cat(4,real(Ehy),imag(Ehy)); % X, Y, Slice, R/I, Echo
        %Ehy = permute(Ehy,[5,4,3,2,1]); % Echo, R/I, Slice, Y, X
        cgsense = cat(4,real(cgsense),imag(cgsense)); % X, Y, Slice, R/I, Echo
        %cgsense = permute(cgsense,[5,4,3,2,1]); % Echo, R/I, Slice, Y, X
        
        for IdxSlice = 1:Nslice
            EHy_SingleSlice = single(squeeze(Ehy(:,:,IdxSlice,:,:)));
            cgsense_SingleSlice = single(squeeze(cgsense(:,:,IdxSlice,:,:)));
            fileName = ['../',Category,'/EHy_R10_2D/sb',num2str(sbj),'/traj_',num2str(Traj_Set),'_sb_',num2str(sbj),'_slice_',num2str(IdxSlice),'.mat']
            save(fileName,'EHy_SingleSlice', 'cgsense_SingleSlice');
            for IdxEcho = 1:Necho
                slice_img = squeeze(EHy_SingleSlice(:,:,1,IdxEcho) + EHy_SingleSlice(:,:,2,IdxEcho)*1j);
                slice_ksp = log(abs(fftshift(fft2(slice_img))));
                slice_ksp = slice_ksp - min(slice_ksp(:));
                slice_ksp = slice_ksp/max(slice_ksp(:));
                slice_img = abs(slice_img);
                slice_img = slice_img/max(slice_img(:));
                PNG2Bsaved = [slice_img,slice_ksp];
                PNG2Bsaved = PNG2Bsaved/max(PNG2Bsaved(:));
                PNGName = ['../',Category,'/EHy_R10_2D/sb',num2str(sbj),'/PNG/traj_',num2str(Traj_Set),'_sb_',num2str(sbj),'_slice_',num2str(IdxSlice),'_Echo_',num2str(IdxEcho),'.png'];
                imwrite(PNG2Bsaved,PNGName);
            end
        end
    end
end

