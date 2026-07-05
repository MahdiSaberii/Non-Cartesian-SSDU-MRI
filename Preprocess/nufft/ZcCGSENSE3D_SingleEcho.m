function  u = ZcCGSENSE3D_SingleEcho(FHy, coil, alpha,RapidOperator,numIter)

%y = FHy/max(abs(FHy(:)));
y = FHy;
[nx,ny,nz,nc] = size(coil);
M  = @(x) applyM(RapidOperator,coil,x) + alpha*x;
x = 0*y(:);
r = y(:);
p = r;
rr = r'*r;

list = [FHy(:,:,40)];
previous = FHy(:,:,40);
difflist = zeros([120,120],'single');
%L2list = [];

for it = 1:numIter
    Ap = M(p);
    a = rr/(p'*Ap);
    x = x + a*p;
    %L2list = [L2list, (r'*r) / (a'*a);];
    %aa=reshape(a*p,nx,ny,nz);
    %figure;imshow(abs(squeeze(aa(:,64,:))),[])
    rnew = r - a*Ap;
    b = (rnew'*rnew)/rr;
    r=rnew;
    rr = r'*r;
    p = r + b*p;
    %disp([num2str(max(abs(col(a)))),';',num2str(max(abs(col(r)))),';',num2str(max(abs(col(b)))),';',num2str(max(abs(col(p)))),';']);

    u_it = reshape(x,nx,ny,nz);
    diff = abs(u_it(:,:,40) - previous)*5;
    list = [list,u_it(:,:,40)];
    difflist = [difflist, diff];
    previous = u_it(:,:,40);
end
%figure;imshow(abs([list;difflist]),[])
u  = reshape(x,nx,ny,nz);
%figure;plot(L2list)
end

%% Derivative evaluation
function y = applyM(M,coil,x)
[nx,ny,slice,nc] = size(coil);
dx = reshape(x,[nx,ny,slice,1]);
ZeroPadded_img = zeros([nx*2,ny*2,slice,nc],'single');
coilimg = repmat(dx,[1,1,1,nc]).*coil;
ZeroPadded_img(1+nx/2:nx/2+nx,1+ny/2:ny/2+ny,:,:)=coilimg;
%size(M)
%size(ZeroPadded_img)
itermediate = repmat(M,[1,1,slice,nc]).*  fft(fft(ZeroPadded_img ,[],1),[],2)/sqrt(2*nx*2*ny);
rapid  = ifft(ifft( itermediate ,[],1),[],2) * sqrt(2*nx*2*ny);
y = (sum(rapid(1+nx/2:nx/2+nx,1+ny/2:ny/2+ny,:,:).*conj(coil),4));
%figure;imshow(abs(y),[]);title('cg intermediate')
y = y(:) ;
end

