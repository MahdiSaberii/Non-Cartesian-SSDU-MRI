training_all =zeros(5,2176);
validation_all = zeros(5, 2176);

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/mask_new_extra_0_4_5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;


%% Calculate NUFFT operator
shift   = [0,0];
Neco = 6;
fprintf('Calculating N operators \n')
%load('/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/RawData/singleShot_sb7.mat')
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

for m = 6: 10
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m-5,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(m, ie,:,:) = M1;
end
end