function  res = NUFFTDCF(k,w,shift,imSize,mask)

%	non uniform 2D fourier transform operator, based on
%	Jeffery Fessler's code.
%
%	Inputs:
%		k - normalized kspace coordinates (complex value between -0.5 to 0.5)
%`		w - density compensation (w=1, means no compensation)
%		shift - shift the image center
%		imSize - the image size
%
%	Outputs:
%		FT = the NUFFT operator
%

om      = [real(k(:)), imag(k(:))]*2*pi;
Nd      = imSize;
Jd      = [6,6];
Kd      = floor(Nd*1.5);
n_shift = Nd/2 + shift;

res.st  = nufft_init(om, Nd, Jd, Kd, n_shift, 'kaiser');

if ~exist('mask','var')
    mask = 'ones';
end

switch mask
    case 'ones'
        res.st.mask = 1;
    case 'circular'
        [xx,yy] = meshgrid(linspace(-0.5,0.5,Kd(1)));
        cutoff = 0.5;
        radius = sqrt(xx.^2 + yy.^2);
        mask = double(radius < cutoff);
        res.st.mask = fftshift(fftshift(mask,1),2);
    case 'smooth'
        [xx,yy] = meshgrid(linspace(-0.5,0.5,Kd(1)));
        cutoff = 0.5;
        radius = sqrt(xx.^2 + yy.^2);
        mask = 0.5 + 1/pi*atan(100*(cutoff-radius)/cutoff);
        res.st.mask = fftshift(fftshift(mask,1),2);
    otherwise
        error('Chosen mask type not supported.');
end

res.adjoint     = 0;
res.imSize      = imSize;
res.dataSize    = size(k);
res.w           = sqrt(w);
res             = class(res,'NUFFTDCF');

