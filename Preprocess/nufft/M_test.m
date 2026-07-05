%% Calculate NUFFT operator
shift   = [0,0];
Neco = 5;
fprintf('Calculating N operators \n')

load('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/dcSingle_2_2_R10E8.mat')
dc10x = permute(dcf, [2,1]);
%%
load ('/home/naxos2-raid7/hongygu/non_cartesian/Paper_data/3_3_single_R10E15.mat')
for ie = 1:Neco
    trj = trj_kxy(:, 1, ie);
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    par.FT{ie} = NUFFTDCF(trj, ones(size(DCF.^2)), shift, [240, 240]);
end

%%
M = zeros(Neco, 240, 240);
for ie = 1:Neco
    DCF = dc10x(:,ie)./ max(max(dc10x(:,ie)));
    W1 = par.FT{ie}' * DCF;
    M1 = fftshift(fft2(ifftshift(W1))) ./ 240;
    M(ie,:,:) = M1;
end
