function Xd = decimateTrans(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,KU] = size(X);
K = KU/numScales;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xd = zeros(N1,N2,K);
% Construct stacked decimation matrices
for k = 1:K
    for j = 1:numScales
        jk = j+(k-1)*numScales;
        jj = 2^(j-1); 
        dj = X(:,1:N2,jk);
        for i = 1:(j-1)
            Ndj = numel(dj);
            dj1 = zeros(1,2*Ndj);
            dj1(:,2:2:2*Ndj) = dj;
            dj = lwpass4(dj1,1);
        end
        Xd(:,:,k) = Xd(:,:,k) + dj(:,1:N2);
        clear dj
    end
end

end

