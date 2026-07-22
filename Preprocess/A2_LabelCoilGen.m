close all; clear all; clc;
addpath("./nufft/")
addpath("./Espirit/")

nEcho               = 5;
nArms               = 6;
TotalNumberofPoints = 3592; 
img_width           = 120;
shift               = [0,0];
data_path           = "./Raw/";
par.prev            = 'single';
par.prec            = "single"; 

if ~exist("./Labels/", 'dir')
        mkdir("./Labels/");
end 

%%
Tukey_filter = load("../Data/Tuke_01.mat").LPfilter;
Tukey_filter = repmat(Tukey_filter, [1 1 1 1 1]);

for subj_idx=11:17
    subj_num = subj_idx;
    
    % Match files like:
    % ref_MES012_4840__MID00023.mat
    Subject_pattern = sprintf("ref_MES%03d_4840__MID*.mat", subj_num);
    file_list       = dir(fullfile(data_path, Subject_pattern));
    name            = fullfile(file_list(1).folder, file_list(1).name);
    
    dir_name_coil  = "./Coils/" + subj_idx;
    dir_name_label = "./Labels/subject_" + subj_idx;
    
    if ~exist(dir_name_coil, 'dir')
        mkdir(dir_name_coil);
    end

    fprintf("Loading subject %d: %s\n", subj_num, name);
    data = load(name); 
    % data shape:  17960,6,24,10,52
    % trj shape :  17960,6

    nCoil     = size(data.data,5); 
    nSlice    = size(data.data,3);
    data_ref  = data.data(:,:,:,1,:);

    load('./DCF/DCF_Full.mat')
    trj_full = kx + 1j * ky;
    for ie = 1:nEcho
        DCF = dcf_full(:,:,ie)./ max(max(dcf_full(:,:,ie)));
        par.FT{ie} = NUFFTDCF(trj_full(:,:,ie), ones(size(DCF.^2)), shift, [img_width, img_width]);
    end

    img_uncombined_full = zeros([[img_width,img_width], nCoil, nEcho, nSlice]);
    data_full_reshaped  = reshape(data_ref, [TotalNumberofPoints, nEcho, nArms, nSlice, nCoil]);

    fprintf('Gridding of k-space data \n')
    for ip = 1:nSlice
        for ie = 1:nEcho
            DCF = dcf_full(:,:,ie)./ max(max(dcf_full(:,:,ie)));
            for ic = 1:nCoil
                fprintf('\t Neco: %d/%d  Ncoil: %d/%d Nslice: %d/%d \n',ie,nEcho,ic,nCoil, ip,nSlice)           
                img_uncombined_full(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_full_reshaped(:,ie,:,ip,ic)), par.prec));          
            end
        end
    end

    img = permute(img_uncombined_full, [1,2,5,3,4]); %[..., ..., nCoil, nEcho, nSlice] -> [..., ..., nSlice, nCoil, nEcho] 
    % Removing corners
    img = img .* Tukey_filter;
    
    CartesianKsp = fftshift(fftshift(fft(fft(fftshift(fftshift(img,1),2),[],1),[],2),1),2);
    
    % Take the first echo for coil senstitivity generation
    CartesianKsp = squeeze(CartesianKsp(:,:,:,:,1));
    [kx_t,ky_t,slices_t,nCh_t]=size(CartesianKsp);
    coil = zeros([img_width,img_width,nSlice,nCoil],'single');
    
    for is=1:nSlice
        % kspace: Multi channel k-space data. Expected dimensions are (sx, sy, nc), where (sx, sy) are volumetric dimensions and (nc) is the channel dimension.
        fprintf('\t ESPIRiT for Nslice: %d/%d \n', is,nSlice)
        temp_ksp       = squeeze(CartesianKsp(:,:,is,:));
        coil(:,:,is,:) = espirit_generator(temp_ksp,6,24,0.02,0.95); 
    end

    %% Saving the Coils and labels

    label_ehy = zeros(img_width,img_width, nEcho, nSlice);
    for i = 1:nSlice
        sprintf("Coils saving... %d", i)
        
        temp_img = squeeze(img(:,:,i,:,:));
        sens_map = squeeze(coil(:,:,i,:));
        label_ehy(:,:,:,i) = squeeze(sum(temp_img .* conj(sens_map), 3));
    
        save_name = dir_name_coil + "/slice_" + num2str(i) + ".mat";
        save(save_name, "sens_map");
    end
    
    for is=1:nSlice
        sprintf("Labels saving... %d", is)
        label_ehy_temp = label_ehy(:,:,:,is);
        save_name_label = dir_name_label + "_slice_" + num2str(is) + ".mat";
        save(save_name_label, "label_ehy_temp");
        
        save_name_label = dir_name_label + "_slice_" + num2str(is) + ".png";
        max_vals = max(abs(label_ehy(:,:,:,is)), [], [1, 2]);  % Size: 1 x 1 x 10
        imwrite(reshape(abs(label_ehy(:,:,:,is)) ./ max_vals, [img_width,img_width*nEcho]), save_name_label);  % Save the 2D slice
    end
end
