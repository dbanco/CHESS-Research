function [ gradTV ] = gradientTV( f, tvBeta, D )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
gradTV = zeros(size(f,2),1);
t = size(f,1) - 1;
psi_prime = 1./sqrt((D*f).^2 + tvBeta^2);
for j = 1:size(f,2)
    L_f = D'*diag(psi_prime(:,j))*D;
    gradTV(j) = L_f(t,:)*f(:,j);
end
    
end

