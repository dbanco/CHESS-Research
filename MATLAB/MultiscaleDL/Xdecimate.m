function Xdecim = Xdecimate(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

[N1,N2,KU,T] = size(X);
K = KU/numScales;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xm = reshape(X,[T,KU,N2*N1]);
Xmdecim = zeros(T,K,N2*N1);

% Construct stacked decimation matrices
A = eye(N2);
for k = 1:K
    for j = 1:numScales
        jj = 2^(j-1);
        i = j+(k-1)*numScales;
        B = zeros(N2);
        B(1:N2/jj,:) = A(jj:jj:N2,:);
        Xmdecim(:,k,:) = Xmdecim(:,k,:)+...
        reshape(squeeze(Xm(:,i,:))*B,[T,N1,N2]);
    end
Xdecim = reshape(Xmdecim,[N1,N2,K,T]);
end
