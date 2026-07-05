function [f,g] = L2reg(x,Psi,lambda)

% L2-norm
% 2D TV

% use TVOP
w = Psi * x;
% use imfilter
% sigma = 2;
% h_filter = fspecial('gaussian', [3 3], sigma);
% w = imfilter(x, h_filter);


% objective
f = lambda * norm(w(:)).^2;

% gradient
g = 2 * w;

end