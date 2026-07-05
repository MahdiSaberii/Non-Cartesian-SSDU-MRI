training_all =zeros(10,2176);
validation_all = zeros(10, 2176);
%%
load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_1.mat')
training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_6.mat')
training_all(6, :) = training;
validation_all(6, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_7.mat')
training_all(7, :) = training;
validation_all(7, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_8.mat')
training_all(8, :) = training;
validation_all(8, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_9.mat')
training_all(9, :) = training;
validation_all(9, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_3_10.mat')
training_all(10, :) = training;
validation_all(10, :) = validation;

%% Calculate NUFFT operator
shift   = [0,0];
Neco = 8;
fprintf('Calculating N operators \n')

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcSingle_2_2_R10E8.mat')
dc10x = permute(dcf, [2,1]);
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_single_R10E8.mat')
for m = 1: 10
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(training_all(m,:)~=1) = 0;
    par.FT{m, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [240, 240]);
end
end

for m = 1: 10
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m,:)~=1) = 0;
    par.FT{m+10, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [240, 240]);
end
end

%%
M = zeros(20, 5, 240, 240);
for m = 1: 10
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(training_all(m,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(m, ie,:,:) = M1;
end
end

for m = 11: 20
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m-10,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(m, ie,:,:) = M1;
end
end