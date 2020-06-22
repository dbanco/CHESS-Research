function DtB = DiffTranB_1D(B)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here
[N,M,T] = size(B);
T = T+1;
DtB = zeros(N,M,T);
DtB(:,:,1) = -B(:,:,1);
DtB(:,:,T) =  B(:,:,T-1);
for t = 2:(T-1)
    DtB(:,:,t) = B(:,:,t-1) - B(:,:,t);
end