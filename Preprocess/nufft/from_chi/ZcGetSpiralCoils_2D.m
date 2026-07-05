function [coil] = ZcGetSpiralCoils(gridding, BallMask)

% if require compression then do it
CartesianKsp = fftshift(fftshift(fft(fft(fftshift(fftshift(gridding,1),2),[],1),[],2),1),2);
[kx,ky,kz,nCh,nEcho]=size(CartesianKsp);

wid = 30/2;

domain = wid*2;
Twin2D = zeros([domain,domain],'single');
tukey1 = tukeywin(domain,1);
tukey1 = tukey1(domain/2+1:domain);
x = linspace(-domain/2, domain/2, domain);
y = linspace(-domain/2, domain/2, domain);
for i=1:1:domain
for j=1:1:domain
  if (round(sqrt(x(i)^2 + y(j)^2)) <= domain/2)
      x(i);
      y(i);
      round(sqrt(x(i)^2+y(j)^2));
  Twin2D(i,j) = tukey1(round(sqrt(x(i)^2+y(j)^2)));
  end
end
end
%{
Twin3D = zeros([wid*2,wid*2,wid*2],'single');
for ii = 1:wid*2
    Twin3D(:,:,ii) = Twin2D.*repmat(squeeze(Twin2D(:,ii)),[1,wid*2]);
end
%}
LPfilter = zeros([kx,ky],'single');
LPfilter(kx/2-wid+1:kx/2+wid,ky/2-wid+1:ky/2+wid) = Twin2D;
coil = zeros([kx,ky,kz,nCh,nEcho],'single');

for sliceIdx = 1: kz
for echo = 1:nEcho
    for ii=1:nCh
        coil(:,:,sliceIdx,ii,echo) = CartesianKsp(:,:,sliceIdx,ii,echo).*LPfilter;
    end
end
end

coil = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(coil,1),2),[],1),[],2),1),2).* repmat(BallMask,[1,1,1,nCh,nEcho]);

coil_sos =  sqrt(sum(abs(coil).^2,4));
coil = coil./(coil_sos+eps*ones(size(coil_sos),'single'));
%coil(find(abs(coil)>2))=0;
coil = single(coil);
end

