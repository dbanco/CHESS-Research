function Xud = upDwnSampleTrans(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,KU] = size(X);
U = 2*numScales+1;
N0 = N2/2^numScales;
K = KU/U;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xud = zeros(N1,N0,K);
Xarray = reshape(X,[N1,N2,K,U]);
% plot(Xarray(:,:,1,2))
Xd = reshape(Xarray(:,:,:,(numScales+1):-1:1),[N1,N2,K*(numScales+1)]);
Xu = reshape(Xarray(:,:,:,(numScales+1):U),[N1,N2,K*(numScales+1)]);
Xtd = decimateTrans(Xd,numScales+1);
Xud = upSampleTrans(Xu,numScales+1);

Xud = Xtd(:,1:N0)+Xud(:,1:N0);
% subtract double counted identity terms
for k = 1:K
    Xud = Xud - X(:,1:N0,numScales+1+U*(k-1));
end


end
