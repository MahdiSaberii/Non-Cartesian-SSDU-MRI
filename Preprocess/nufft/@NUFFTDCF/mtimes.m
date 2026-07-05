function ress = mtimes(a,bb)

for n = 1:size(bb,3)
    b = bb(:,:,n);
    
    if a.adjoint
        
        b   = b(:).*a.w(:);
        res = nufft_adj(b, a.st) / sqrt(prod(a.imSize));
        res = reshape(res, a.imSize(1), a.imSize(2));
        
        % Account for undersampling
        % res = res .* size(a.w,1) * pi / 2 / size(a.w,2);
        
    else
        
        b   = reshape(b, a.imSize(1),a.imSize(2));
        res = nufft(b, a.st) / sqrt(prod(a.imSize));
        res = res .* a.w(:);
        res = reshape(res,a.dataSize(1),a.dataSize(2));
        
    end
    ress(:,:,n) = res;
end
