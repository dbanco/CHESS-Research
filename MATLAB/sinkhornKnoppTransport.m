function [Wd, r, c, T] = sinkhornKnoppTransport(D,lam,r,c)
%sinkhornKnoppTransport

iter = 1;
Ind = r > 0;
r = r(Ind);
D = D(Ind,:);
K = exp(-lam*D);

x = ones(size(r));
x_old = x/2;

while 1
    e = c./(K'*(1./x));
    x = diag(1./r)*K*e;

    check = e.*K'*(1./x) - c;
    if sum(abs(check ))<1e-8
        break
    end
    x_old = x;
end

u = 1./x; 
v = c.*(1./(K'*u));
T = diag(u)*K*diag(v);
Wd = sum(u.*((K.*D)*v));

end
