function  res = MCNUFFT(k,w,b1,mask)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler and the single-coil NUFFT
% operator from Miki Lustig
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng & Ricardo Otazo, NYU, 2012

Nd = size(b1(:,:,1));
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;

%tbenkert
if ~exist('mask','var')
    mask = 'ones';
end

switch mask
    case 'ones'
        mask_used = 1;
    case 'circular'
        [xx,yy] = meshgrid(linspace(-0.5,0.5,Kd(1)));
        cutoff = 0.5;
        radius = sqrt(xx.^2 + yy.^2);
        mask = double(radius < cutoff);
        mask_used = fftshift(fftshift(mask,1),2);
    case 'smooth'
        [xx,yy] = meshgrid(linspace(-0.5,0.5,Kd(1)));
        cutoff = 0.5;
        radius = sqrt(xx.^2 + yy.^2);
        mask = 0.5 + 1/pi*atan(100*(cutoff-radius)/cutoff);
        mask_used = fftshift(fftshift(mask,1),2);
    otherwise
        error('Chosen mask type not supported.');
end

if isempty(w)
    calcDCF = 1;
end

for tt=1:size(k,3),
    kk=k(:,:,tt);
    om = [real(kk(:)), imag(kk(:))]*2*pi;
    res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
    res.st{tt}.mask = mask_used;
    
    if calcDCF
        w_tmp = my_dcf(res.st{1});
        w(:,:,tt) = reshape(w_tmp,[size(k,1), size(k,2)]);
        % figure(303), plot(reshape(w,[size(k,1), size(k,2)]))
    end
    
end

res.adjoint = 0;
res.imSize = size(b1(:,:,1));
res.dataSize = size(k);
res.w = sqrt(w);
res.b1 = b1;
res = class(res,'MCNUFFT');

