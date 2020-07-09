function DtPtB = DiffTranB_1D(B)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here
[N,K,T] = size(B);
T = T+1;
DtPtB = zeros(N,K,T);
DtPtB(:,:,1) = -B(:,:,1);
DtPtB(:,:,T) =  B(:,:,T-1);
for t = 2:(T-1)
    DtPtB(:,:,t) = B(:,:,t-1) - B(:,:,t);
end