function [f,g] = L1reg(x,Psi,lambda,l1Smooth)

% L1-norm
% 2D TV

w = Psi * x;

% objective
f = lambda * sum((conj(w(:)).*w(:) + l1Smooth).^(1/2));

% gradient
g = lambda * (Psi'*(real(w).*(w.*conj(w)+l1Smooth).^(-0.5)));

end