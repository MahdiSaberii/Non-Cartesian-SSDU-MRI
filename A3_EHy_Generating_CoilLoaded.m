close all; clear all; clc;

N_Masks             = 7;
nEcho               = 5;
nArms               = 6;
Acc_R               = nArms; 
TotalNumberofPoints = 3592;  
img_width           = 120;
shift               = [0,0];
par.prev            = 'single';
par.prec            = "single";
data_path           = "./Raw/Raw_Lite/";
saving_path         = sprintf("./Ehy_%d", Acc_R);
time_idxs           = [20, 55, 90, 125]; % fMRI time frames to take
Train_Subjects      = [11,12,13,14];     % Subject numbers to take

if ~exist(saving_path, 'dir')
        mkdir(saving_path);
end

%% Loading Multi-Mask SSDU Masks
training_all   = zeros(N_Masks,TotalNumberofPoints);
validation_all = zeros(N_Masks,TotalNumberofPoints);
for i=1:N_Masks
    name = "./Masks/mask_"+num2str(i)+".mat";
    load(name)
    training_all(i, :)  = training;
    validation_all(i,:) = validation;
end


%%
for subj_idx=1:length(Train_Subjects)
    subj_num       = Train_Subjects(subj_idx);
    % Match files like:
    % MES012_4840__MID00023_4TF.mat
    Subject_pattern = sprintf("MES%03d_4840__MID*_4TF.mat", subj_num);
    file_list       = dir(fullfile(data_path, Subject_pattern));
    name            = fullfile(file_list(1).folder, file_list(1).name);

    data   = load(name); % shape: 17960, 24, 4, 52 --> (Spikes*Echo, Slice, TimeFrames, Coils)
    nSlice = size(data.data, 2);
    nCoil  = size(data.data, 4); 
%% Loading coil sensitivity maps
    coil = zeros([img_width,img_width,nSlice,nCoil],'single');
    for is=1:nSlice
        % kspace: Multi channel k-space data. Expected dimensions are (sx, sy, nc), where (sx, sy) are volumetric dimensions and (nc) is the channel dimension.
        fprintf('\t ESPIRiT for Nslice: %d/%d \n', is,nSlice)

        filename       = sprintf("./Coils/%d/slice_%d.mat", subj_num, is);
        data_temp      = load(filename);  % shape: 120*120*52
        coil_temp      = data_temp.sens_map;
        coil(:,:,is,:) = coil_temp; 
    end

%% Loading the DCF function
    dcf_single = load('./DCF/DCF_Single.mat');
    trj_single = dcf_single.kx + 1j * dcf_single.ky;
    for ie = 1:nEcho
        DCF = dcf_single.dcf_single(:,ie)./ max(max(dcf_single.dcf_single(:,ie)));
        par.FT{ie} = NUFFTDCF(trj_single(:,ie), ones(size(DCF.^2)), shift, [img_width, img_width]);
    end
%% Singleshot Preprocessing
    data_single_huge = data.data; 
    img_uncombined   = zeros([[img_width,img_width], nCoil, nEcho, nSlice]);
    
    for idx_time = 1:length(time_idxs)
        data_single = data_single_huge(:,:,idx_time, :);
        data_single_reshaped = reshape(data_single, [TotalNumberofPoints, nEcho, nSlice, nCoil]);

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
        img_uncombined            = permute(img_uncombined, [1,2,5,3,4]);% -> 120,120,24,34,10
        img_uncombined_train      = zeros([[img_width,img_width], nCoil,nEcho,nSlice]);
        img_uncombined_validation = zeros([[img_width,img_width], nCoil, nEcho, nSlice]);
        training_data             = zeros(N_Masks, TotalNumberofPoints, nEcho, nSlice, nCoil);
        validation_data           = zeros(N_Masks,TotalNumberofPoints, nEcho, nSlice, nCoil);

        for m = 1: N_Masks
            training_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
            training_data(m, training_all(m,:)~=1,:,:,:) = 0;
            validation_data(m, :,:,:,:) = data_single_reshaped(:, :, :,:,:);
            validation_data(m, training_all(m,:)==1,:,:,:) = 0;
        end
%%
        fprintf('Gridding of k-space data \n')
        theta_ehy  = zeros(img_width,img_width,nEcho, nSlice, N_Masks);
        lambda_ehy = zeros(img_width,img_width,nCoil,nEcho, nSlice, N_Masks);

        for m=1:N_Masks
        for is = 1:nSlice
            for ie = 1:nEcho
                DCF_train = transpose(dcf_single.dcf_single(:,ie)./ max(max(dcf_single.dcf_single(:,ie)))); 
                DCF_val = transpose(dcf_single.dcf_single(:,ie)./ max(max(dcf_single.dcf_single(:,ie))));
                DCF_train(training_all(m,:)~=1) = 0;
                DCF_val(validation_all(m,:)~=1) = 0;
                trj_train = trj_single(:,ie);
                trj_val = trj_single(:,ie);
                trj_train(training_all(m,:)~=1, :)=0;
                trj_val(validation_all(m,:)~=1, :)=0;
                FT_train = NUFFTDCF(trj_train,ones(size(DCF.^2)), shift, [120, 120]);
                FT_val = NUFFTDCF(trj_val, ones(size(DCF.^2)), shift, [120, 120]); 
                for ic = 1:nCoil
                    fprintf('\t Echo#: %d/%d  Coil#: %d/%d Slice#: %d/%d  Mask#: %d/%d  \n',ie,nEcho,ic,nCoil, is, nSlice, m, N_Masks)           
                    img_uncombined_train(:,:,ic,ie, is) = FT_train' * (DCF_train.*cast(squeeze(training_data(m,:,ie,is,ic)), par.prec));          
                    img_uncombined_validation(:,:,ic,ie,is) = FT_val' * (DCF_val.*cast(squeeze(validation_data(m,:,ie,is,ic)), par.prec));
                end
            end

            theta_ehy(:,:,:,is, m)   = sum(img_uncombined_train(:,:,:,:,is).*conj(squeeze(coil(:,:,is,:))), 3);
            lambda_ehy(:,:,:,:,is,m) = img_uncombined_validation(:,:,:,:,is);
        end
        end

%%
        dir_name = sprintf("./Ehy_%d/subject_%d",Acc_R, subj_num);
        
        for is=1:nSlice
            fprintf('\t Saving Slice: %d/%d, Time: %d \n',is,nSlice, time_idxs(idx_time))
            theta_ehy_6  = squeeze(theta_ehy(:,:,:,is,:));
            lambda_ehy_6 = squeeze(lambda_ehy(:,:,:,:,is,:));
            save_name = dir_name + "_slice_" + num2str(is)+"_t_"+num2str(time_idxs(idx_time))+ ".mat";
            save(save_name, "theta_ehy_6", "lambda_ehy_6");
        
            save_name_ehy   = dir_name + "_slice_" + num2str(is) +"_t_"+num2str(time_idxs(idx_time))+ ".png";
            theta_ehy_6_png = squeeze(theta_ehy(:,:,:,is,1)); % Taking the first mask data, for visualization
            max_vals        = max(abs(theta_ehy_6_png), [], [1, 2]);  % Size: 1 x 1 x 10
            
            ehy_combined = reshape(theta_ehy_6_png./max_vals,[img_width,img_width*nEcho]);
            imwrite(abs(ehy_combined), save_name_ehy);  % Save the 2D slice
        end
    end
end