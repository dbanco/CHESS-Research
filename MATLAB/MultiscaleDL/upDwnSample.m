function Dud = upDwnSample(D,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,K] = size(D);
U = 2*numScales+1;
nSp1 = numScales+1;
Dud = zeros(N1,N2*(2^numScales),U,K);
for k = 1:K
    Du = upSample(D(:,:,k),nSp1);
    Dd = decimate(D(:,:,k),nSp1);

    Dud(:,1:N2,1:numScales,k) = Dd(:,:,nSp1:-1:2);
    Dud(:,1:N2,nSp1,k) = Dd(:,:,1);
    Dud(:,:,(nSp1+1):U,k) = Du(:,:,2:nSp1);
end
Dud = reshape(Dud,[N1,N2*(2^numScales),K*U]);




