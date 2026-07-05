
function wi = my_dcf(G)
niter   = 6;
thresh  = 0.02;

m       = G.M;
N       = G.Nd;
k       = G.om/(2*pi);
P       = G.p;

wi      = ones(m,1);
wi_filt = ones(m,1);

goal    = inf;
iter    = 0;
saver   = zeros(niter,1);

while max(abs(goal-1)) > thresh
    iter = iter + 1;
    
    
    if isfield(G,'interp_table')
        goal = feval( G.interp_table, G, feval(G.interp_table_adj, G, wi) );
    else
        goal = P*(P'*wi);
    end
    
    wi = (wi.*wi_filt) ./ abs(goal);
    
    if iter > niter
        warning 'iteration stuck?'
        break
    end
    saver(iter) = max(abs(goal-1));
end
%

% filtering by fxbreuer

r      = sqrt(sum(k.^2,2));
rMax   = max(r(:));
idx    = find(r > 0.9*rMax);

wi_filt = medfilt1(wi,9);
wi(idx) = wi_filt(idx);

% scaling by fxbreuer
wi = wi./sum(abs(wi(:)));
wi = wi * (rMax).^2*pi;
wi = prod(N)*wi;

printm('pipe ended at iteration %d with %g', iter, max(abs(goal-1)))

% function wi = mri_dcf_samsonov(G)
% 
% niter   = 20;
% thresh  = 0.02;
% 
% m       = G.M;
% N       = G.Nd;
% k       = G.om/(2*pi);
% P       = G.p;
% 
% wi      = ones(m,1);
% wi_filt = ones(m,1);
% 
% goal    = inf;
% iter    = 0;
% saver   = zeros(niter,1);
% 
% R      = sqrt(sum(k.^2,2));
% rMax   = max(R(:));
% idx    = find(R > 0.9*rMax);
% 
% iter = iter + 1;
% 
% if isfield(G,'interp_table')
%     tmp = feval( G.interp_table, G, feval(G.interp_table_adj, G, wi) );
% else
%     tmp = P * (P' * wi);
% end
% %wi = (1 ./ tmp);
% wi = real(1 ./ tmp);
% wi0=wi;
% 
% for iter=1:niter
%     
%     
%     if isfield(G,'interp_table')
%         tmp = feval( G.interp_table, G, feval(G.interp_table_adj, G, wi) );
%     else
%         tmp = P * (P' * wi);
%     end
%     
%     r = wi0.*(ones(m,1)-tmp);
%     if isfield(G,'interp_table')
%         tmp = feval( G.interp_table, G, feval(G.interp_table_adj, G, r) );
%     else
%         tmp = P * (P' * r);
%     end
%     alpha = (r'*r)./(r'*tmp);
%     
%     
%     wi(wi<0)=0;
%     wi = real(wi + alpha*r);
%     %       wi = wi + alpha*r;
%     figure(100),plot(real(wi))
% end
% 
% 
% % filtering by fxbreuer
% 
% r      = sqrt(sum(k.^2,2));
% rMax   = max(r(:));
% idx    = find(r > 0.9*rMax);
% 
% 
% wi_filt = medfilt1(wi,9);
% %wi(idx) = wi_filt(idx);
% 
% % scaling by fxbreuer
% wi = wi./sum(abs(wi(:)));
% wi = wi * (rMax).^2*pi;
% wi = prod(N)*wi;
% 
% printm('samsonov ended at iteration %d with %g', iter, max(abs(goal-1)))

    

% function wi = mri_dcf_bydder(G)
% 
% niter   = 6;
% thresh  = 0.02;
% 
% m       = G.M;
% N       = G.Nd;
% k       = G.om/(2*pi);
% P       = G.p;
% 
% wi      = ones(m,1);
% wi_filt = ones(m,1);
% 
% goal    = inf;
% iter    = 0;
% saver   = zeros(niter,1);
% beta    = 1e2;
% 
% n=prod(G.Kd);
% 
% ones_filt_cart  = ones(n,1);
% ones_filt       = ones(m,1);
% 
% saver = zeros(niter,1);
% 
% 
% if isfield(G,'interp_table')
%     tmp = feval( G.interp_table, G, feval(G.interp_table_adj, G, wi) );
% else
%     tmp = P*(P'*wi);
% end
% 
% wi = (1 ./ tmp);
% g = 1./wi + beta.^2;
% 
% if isfield(G,'interp_table')
%     r = feval( G.interp_table, G, ones_filt_cart-feval(G.interp_table_adj, G, wi) );
%     r= ones_filt.*r;
% else
%     r=P*(ones_filt_cart-P'*wi);
% end
% s = r./g;
% new = r'*s;
% 
% %figure(100),plot(abs(wi(1:2000)))
% for iter=1:niter
%     
%     if isfield(G,'interp_table')
%         q = ones_filt.*(feval( G.interp_table, G, feval(G.interp_table_adj, G, s) )) + beta.^2.*s;
%     else
%         q = P*(P'*s)+beta.^2.*s;
%     end
%     alpha = new/(s'*q);
%     wi=wi+ones_filt.*(alpha*s);
%     if(sum(wi<0)>0)
%         sprintf('negative values in the weights... may want to try increasing beta')
%         wi(wi<0)=0;
%         %beta=beta*1.1
%         %wi(wi<0)=0;
%     end
%     % figure,plot(real(wi(1:2000)))
%     r=r-alpha*q;
%     q=r./g;
%     old=new;
%     new=r'*q
%     
%     s=q+(new/old)*s;
% end
% 
% % filtering by fxbreuer
% 
% r      = sqrt(sum(k.^2,2));
% rMax   = max(r(:));
% idx    = find(r > 0.9*rMax);
% 
% wi_filt = medfilt1(wi,9);
% wi(idx) = wi_filt(idx);
% 
% % scaling by fxbreuer
% wi = wi./sum(abs(wi(:)));
% wi = wi * (rMax).^2*pi;
% wi = prod(N)*wi;
% 
% 
% 
% wi = abs(wi);
% 
% printm('bydder ended at iteration %d', iter)
% 
