num = 1488;
Neco = 12;

training_all =zeros(10,num);
validation_all = zeros(10, num);
%%
load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_1.mat')

training_all(1, :) = training;
validation_all(1, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_2.mat')
training_all(2, :) = training;
validation_all(2, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_3.mat')
training_all(3, :) = training;
validation_all(3, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_4.mat')
training_all(4, :) = training;
validation_all(4, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_5.mat')
training_all(5, :) = training;
validation_all(5, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_6.mat')
training_all(6, :) = training;
validation_all(6, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_7.mat')
training_all(7, :) = training;
validation_all(7, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_8.mat')
training_all(8, :) = training;
validation_all(8, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_9.mat')
training_all(9, :) = training;
validation_all(9, :) = validation;

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/masks/mask_paper_5_10.mat')
training_all(10, :) = training;
validation_all(10, :) = validation;

%% Calculate NUFFT operator
shift   = [0,0];
fprintf('Calculating N operators \n')

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcSingle_3_3_R8E12.mat')
dc10x = permute(dcf, [2,1]);
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_single_R8E12.mat')
for m = 1: 10
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(training_all(m,:)~=1) = 0;
    par.FT{m, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [160, 160]);
end
end

for m = 1: 10
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m,:)~=1) = 0;
    par.FT{m+10, ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [160, 160]);
end
end

%%
M = zeros(20, Neco, 160, 160);
for m = 1: 10
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(training_all(m,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 160;
    M(m, ie,:,:) = M1;
end
end

for m = 11: 20
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    DCF(validation_all(m-10,:)~=1) = 0;
    W1 = par.FT{m, ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 160;
    M(m, ie,:,:) = M1;
end
end