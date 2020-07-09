 function PtDtR = PhiTranDiffTran_1D(R,N)
%conjGradResidual Compute residual for conjugate gradient that includes 
% difference matrix
%   Detailed explanation goes here
[K,T] = size(R);
T = T+1;
DtR = zeros(K,T);
DtR(:,1) = -R(:,1);
DtR(:,T) =  R(:,T-1);
for t = 2:(T-1)
    DtR(:,t) = R(:,t-1) - R(:,t);
end

% PtDtR = zeros(N,K,T);
% for i = 1:T
%     vec = DtR(:,i)';
%     bigVec = repmat(vec,[N,1]);
%     PtDtR(:,:,i) = bigVec;
% end
DtR = reshape(DtR,1,size(DtR,1),size(DtR,2));
PtDtR = repmat(DtR, [N,1,1]);
