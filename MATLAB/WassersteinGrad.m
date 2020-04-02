function [ gradW, Wd, T ] = WassersteinGrad(r, c, lam, D)
%WassersteinGrad Summary of this function goes here
%   Detailed explanation goes here

iter = 1;

Ind = r > 0;
r = r(Ind);
r = r./sum(r(:));

D = D(Ind,:);
K = exp(-D*lam-1);

x = ones(size(r))/2;
k = 0;

while (k < 1000) && (numel(r) > 0)

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
gradW = zeros(size(c));
gradW(Ind) = (log(u)-log(sum(u)./numel(c)))/lam;

