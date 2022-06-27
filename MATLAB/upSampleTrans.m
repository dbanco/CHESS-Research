function Xdo = upSampleTrans(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

[N1,N2,KU] = size(X);
K = KU/numScales;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xd = zeros(N1,N2,KU/numScales);
% Construct stacked decimation matrices
for k = 1:K
    for j = 1:numScales
        jk = j+(k-1)*numScales;
        dj = X(:,1:N2,jk);
        for i = 1:(j-1)
            dj1 = lwpass4(dj,1);
            dj = dj1(2:2:end);
        end
        Ndj = numel(dj);
        Xd(:,1:Ndj,k) = Xd(:,1:Ndj,k) + dj;
        clear dj
    end
end
Xdo = Xd(:,1:N2/2^(numScales-1),:);
end

