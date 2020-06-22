function Diff = DiffX_1D(X)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here

[N,M,T] = size(X);
Diff = zeros(N,M,T-1);
for t = 1:(T-1)
    Diff(:,:,t) = X(:,:,t+1)-X(:,:,t);
end

