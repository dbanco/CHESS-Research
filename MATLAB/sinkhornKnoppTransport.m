function [Wd, r, c, T] = sinkhornKnoppTransport(r, c, lam, D)
%sinkhornKnoppTransport

iter = 1;

% Remove 0 probability points from each distribution
Indr = r > 0;
r = r(Indr);
r = r./sum(r(:));

Indc = c > 0;
c = c(Indc);
c = c./sum(c(:));

D = D(Indr,Indc);
K = exp(-lam*D-1);

x = ones(size(r))/2;

while 1 && (numel(r) > 0)
    e = c./(K'*(1./x));
    x = diag(1./r)*K*e;

    check = e.*K'*(1./x) - c;
    if sum(abs(check ))<1e-8
        break
    end
end

u = 1./x; 
v = c.*(1./(K'*u));
T = diag(u)*K*diag(v);
Wd = sum(u.*((K.*D)*v));

end
