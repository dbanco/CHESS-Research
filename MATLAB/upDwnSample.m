function Xud = upDwnSample(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,K,T] = size(X);
U = 2*numScales+1;
nSp1 = numScales+1;
Xud = zeros(N1,N2*(2^numScales),K*U,T);
Xu = upSample(X,nSp1);
Xd = decimate(X,nSp1);

Xud(:,1:N2,1:numScales,:) = Xd(:,:,nSp1:-1:2,:);
Xud(:,1:N2,nSp1,:) = Xd(:,:,1,:);
Xud(:,:,(nSp1+1):U,:) = Xu(:,:,2:nSp1,:);




