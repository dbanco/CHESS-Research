function [ gradW ] = WassersteinGrad( h, y, lam, D )
%WassersteinGrad Summary of this function goes here
%   Detailed explanation goes here

Ind = h > 0;
h = h(Ind);
D = D(Ind,:);
K = exp(-lam*D);
u = ones(size(h));
u_old = u/2;


while 1
    u = h./(K*(y./(K'*u)));
    
    if norm(u-u_old) < 1e-8
       break 
    end
    u_old = u;
end
gradW = zeros(size(y));
gradW(Ind) = lam*log(u)-sum(lam*log(sum(u))./numel(h),2);

