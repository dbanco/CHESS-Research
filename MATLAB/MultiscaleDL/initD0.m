function D0 = initD0(K,M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
D0 = zeros(1,M,K); 
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
center = (M+1)/2;
widths = floor(linspace(0,0.5*M,K+1));
for k = 1:K   
    D0(1,center-widths(k+1):center+widths(k+1),k) = 1;
end

D0 = Pnrm(D0);