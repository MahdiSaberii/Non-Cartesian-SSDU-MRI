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

Neco = par.Neco;
Ncoil = par.NcoilFinal;
Npar = par.Npar;

%% Calculate NUFFT operator
shift   = [0,0];
fprintf('Calculating NUFFT operators \n')
for ie = 1:Neco
    par.FT{ie} = NUFFTDCF(traj(:,:,ie,1), dcf(:,:,1), shift, par.fftsize, par.mask);
end

%% Gridding reconstruction
under_mkdata = zeros(512, 600, 4, 6, 44);
under_mkdata(:, 1:4:end, :, :, :) = mkdata(:, 1:4:end, :, :, :);

img_uncombined = zeros([par.fftsize, Ncoil, Neco, Npar]);
fprintf('Gridding of k-space data \n')
for ip = 1:Npar
    for ie = 1:Neco
        for ic = 1:Ncoil
            if par.verbose
                fprintf('\t Neco: %d/%d  Ncoil: %d/%d \n',ie,Neco,ic,Ncoil)
            end            
            img_uncombined(:,:,ic,ie,ip) = par.FT{ie}' * cast(squeeze(under_mkdata(:,:,ic,ie,ip)), par.prec);          
        end
    end
end
% img = img_uncombined(:,:,1,1,1);
% 
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]
% 
% 
% figure, show3d(abs(im(:,:,:,1)), [0 max(abs(im(:)))/10]);
% 
% figure, imshow(abs(im(:,:,1,1)), [])



%%
% img_gridding_sos = squeeze(sqrt(sum(abs(img_uncombined).^2,3)));
% im = permute(squeeze(img_gridding_sos),[1 2 4 3]); % im [3 dims, echoes]

