close all; clear all; clc;

N_Masks             = 1;
nEcho               = 5;
nArms               = 6;
Acc_R               = nArms; 
TotalNumberofPoints = 3592; 
img_width           = 120;
shift               = [0,0];
time_idxs           = 1:136;
Test_Subjects       = [15,16,17];

par.prev    = 'single';
par.prec    = 'single'; 
data_path   = "./Raw/";
saving_path = sprintf("./Ehy_%d_fMRI", Acc_R);
if ~exist(saving_path, 'dir')
        mkdir(saving_path);
end
%%
for subj_idx=1:length(Train_Subjects)
    subj_num       = Train_Subjects(subj_idx);
    
    % Match files like:
    % MES012_4840__MID00023.mat
    Subject_pattern = sprintf("MES%03d_4840__MID*.mat", subj_num);
    file_list       = dir(fullfile(data_path, Subject_pattern));
    name            = fullfile(file_list(1).folder, file_list(1).name);
    
    data             = load(name); % shape: 17960, 24, 136, 52 --> (TotalNumberofPoints*Echo, Slice, TimeFrames, Coils)
    data_single_huge = data.data; 
    nSlice           = size(data.data, 2);
    nCoil            = size(data.data, 4);
    
    coil = zeros([img_width,img_width,nSlice,nCoil],'single');
    for is=1:nSlice
        % kspace: Multi channel k-space data. Expected dimensions are (sx, sy, nc), where (sx, sy) are volumetric dimensions and (nc) is the channel dimension.
        fprintf('\t ESPIRiT for Nslice: %d/%d \n', is,nSlice)

        filename = sprintf("./Coils/%d/slice_%d.mat", subj_num, is);
        data_temp = load(filename);  % Load the .mat file
        coil_temp = data_temp.sens_map;
        coil(:,:,is,:) = coil_temp ; 
    end

%%
    dcf_single = load('../Data/DCF_Single.mat');
    trj_single = dcf_single.kx + 1j * dcf_single.ky;
    for ie = 1:nEcho
        DCF = dcf_single.dcf_single(:,ie)./ max(max(dcf_single.dcf_single(:,ie)));
        par.FT{ie} = NUFFTDCF(trj_single(:,ie), ones(size(DCF.^2)), shift, [img_width, img_width]);
    end

    img_uncombined       = zeros([[img_width,img_width], nCoil, nEcho, nSlice]);
    for idx_time = time_idxs
        data_single = data_single_huge(:,:,idx_time, :);
        data_single_reshaped = reshape(data_single, [totalpoints, nEcho, nSlice, nCoil]);

        fprintf('Gridding of k-space data \n')
        for ip = 1:nSlice
            for ie = 1:nEcho
                DCF = dcf_single.dcf_single(:,ie)./ max(max(dcf_single.dcf_single(:,ie)));
                for ic = 1:nCoil
                    fprintf('\t Echo#: %d/%d  Coil#: %d/%d Slice#: %d/%d \n',ie,nEcho,ic,nCoil, ip,nSlice)          
                    img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * (DCF.*cast(squeeze(data_single_reshaped(:,ie,ip,ic)), par.prec));          
                    % 120,120,34,10,24
                end
            end
        end
        img_uncombined = permute(img_uncombined, [1,2,5,3,4]);% -> 120,120,24,34,10
        img_uncombined_train = zeros([[img_width,img_width], nCoil,nEcho,nSlice]);
        training_data        = zeros(N_Masks, totalpoints, nEcho, nSlice, nCoil);

        for m = 1: N_Masks
            training_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
        end
%%
        fprintf('Gridding of k-space data \n')
        theta_ehy  = zeros(img_width,img_width,nEcho, nSlice, N_Masks);

        for m=1:N_Masks
        for is = 1:nSlice
            for ie = 1:nEcho
                DCF_train = transpose(dcf_single.dcf_single(:,ie)./ max(max(dcf_single.dcf_single(:,ie)))); \
                trj_train = trj_single(:,ie);
                FT_train = NUFFTDCF(trj_train,ones(size(DCF.^2)), shift, [120, 120]);
                
                for ic = 1:nCoil
                    fprintf('\t Echo#: %d/%d  Coil#: %d/%d Slice#: %d/%d  Mask#: %d/%d  \n',ie,nEcho,ic,nCoil, is, nSlice, m, N_Masks)           
                    img_uncombined_train(:,:,ic,ie, is) = FT_train' * (DCF_train.*cast(squeeze(training_data(m,:,ie,is,ic)), par.prec));          
                end
            end
            theta_ehy(:,:,:,is, m)   = sum(img_uncombined_train(:,:,:,:,is).*conj(squeeze(coil(:,:,is,:))), 3);
            
        end
        end

%%
        dir_name = "./Ehy_6_fMRI/subject_"  + subj_num;
        for is=1:nSlice
            fprintf('\t Saving Slice: %d/%d, Time: %d \n',is,nSlice, idx_time)
            theta_ehy_6  = theta_ehy(:,:,:,is,:);
            save_name = dir_name + "_slice_" + num2str(is)+"_t_"+num2str(idx_time)+ ".mat";
            save(save_name, "theta_ehy_6");
        
            save_name_ehy   = dir_name + "_slice_" + num2str(is) +"_t_"+num2str(idx_time)+ ".png";
            theta_ehy_6_png = theta_ehy(:,:,:,is,1); % Taking the first mask data, for visualization
            max_vals        = max(abs(theta_ehy_6_png), [], [1, 2]);  % Size: 1 x 1 x 10
            
            ehy_combined = reshape(theta_ehy_6_png./max_vals,[img_width,img_width*nEcho]);
            imwrite(abs(ehy_combined), save_name_ehy);  % Save the 2D slice
        end
    end
end
