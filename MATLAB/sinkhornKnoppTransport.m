function [Wd, r, c, T] = sinkhornKnoppTransport(r, c, lam, D)
%sinkhornKnoppTransport

% Remove 0 probability points from each distribution
Ind = r > 0;
r = r(Ind);
r = r./sum(r(:));
c = c./sum(c(:));

D = D(Ind,:);
K = exp(-D*lam-1);

x = ones(size(r))/2;

k = 0;
while numel(r) > (k < 1000)

    e = c./(K'*(1./x));
    x = diag(1./r)*K*e;

    check = e.*K'*(1./x) - c;
    if sum(abs(check ))<1e-8
        break
    end
    k = k + 1;

end

u = 1./x; 
v = c.*(1./(K'*u));
T = diag(u)*K*diag(v);
Wd = sum(u.*((K.*D)*v));

end
