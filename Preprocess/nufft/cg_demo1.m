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

a = par.FT{1}' * cast(squeeze(mkdata(:,:,2,1,1)), par.prec);
b = par.FT{1} * a;
c = par.FT{1}' * b;
figure, imshow(abs(c), [])

%%
