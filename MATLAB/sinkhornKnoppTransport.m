function [Wd, r, c, T] = sinkhornKnoppTransport(M,lam,r,c)
%sinkhornKnoppTransport

iter = 1;
Ind = (r > 0);
r = r(Ind);
M = M(Ind,:);
K = exp(-lam*M);

x = ones(numel(r),1)/length(r);
x_old = ones(size(x));

% subsequent iterations include test
while 1
    iter = iter + 1;
    x = diag(1./r)*K*(c.*(1./(K'*(1./x))));
    if norm(x_old - x) <1e-8
        break
    end
    x_old = x;
end

u = 1./x; 
v = c.*(1./(K'*u));
T = diag(u)*K*diag(v);
Wd = sum(u.*((K.*M)*v));

end
