function [ gradW ] = WassersteinGrad( r, c, lam, D )
%WassersteinGrad Summary of this function goes here
%   Detailed explanation goes here

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
gradW = zeros(size(c));
gradW(Ind) = lam*log(u)-sum(lam*log(sum(u))./numel(c),2);

