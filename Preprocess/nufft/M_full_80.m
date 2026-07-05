%% Calculate NUFFT operator
shift   = [0,0];
Neco = 10;
fprintf('Calculating N operators \n')

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcFull_3_3_R6E10.mat')
%%
dc_full = permute(dcf, [3,2,1]);
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_full_R6E10.mat')

%%
for ie = 1:Neco
    trj = trj_kxy(:, :, 1, ie);
    DCF = dc_full(:,:,ie)./ max(max(dc_full(:,:,ie)));
    par.FT{ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [160, 160]);
end

%%
M = zeros(Neco, 160, 160);
for ie = 1:Neco
    DCF = dc_full(:,:,ie)./ max(max(dc_full(:,:,ie)));
    W1 = par.FT{ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 160;
    M(ie,:,:) = M1;
end

%%
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_full_R6E10.mat')
a = squeeze(trj_kxy(:,1,1,:));

%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_single_R6E10.mat')
b = squeeze(trj_kxy(:,1,:));


%%
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_full_R6E5.mat')
a = squeeze(trj_kxy(:,1,1,:));

%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/2_2_single_R6E5.mat')
b = squeeze(trj_kxy(:,1,:));


