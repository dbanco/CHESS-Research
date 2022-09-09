function Xudo = upDwnSampleTrans(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,KU] = size(X);
U = 2*numScales+1;
N0 = N2/2^numScales;
K = KU/U;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xud = zeros(N1,N2,K);
Xarray = reshape(X,[N1,N2,U,K]);
% plot(Xarray(:,:,1,2))
for k = 1:K
    Xd = reshape(Xarray(:,:,(numScales+1):-1:1,k),[N1,N2,(numScales+1)]);
    Xu = reshape(Xarray(:,:,(numScales+1):U,k),[N1,N2,(numScales+1)]);
    Xtd = decimateTrans(Xd,numScales+1);
    Xtu = upSampleTrans(Xu,numScales+1);

    Xud(:,1:N0,k) = Xud(:,1:N0,k) + Xtd(:,1:N0,:) + Xtu...
        - Xarray(:,1:N0,numScales+1,k);
    % subtract double counted identity terms
     Xudo = Xud(:,1:N0,:);
end


end
