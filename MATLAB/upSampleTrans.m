function Xd = upSampleTrans(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
[N1,N2,KU] = size(X);
K = KU/numScales;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xd = zeros(N1,N2/2^(numScales-1),KU/numScales);
% Construct stacked decimation matrices
for k = 1:K
    for j = 1:numScales
        jk = j+(k-1)*numScales;
        dj = X(:,1:N2/2^(numScales-j),jk);
        for i = 1:(j-1)
            dj1 = lwpass4(dj,1);
            dj = dj1(2:2:end);
        end
        Xd(:,:,k) = Xd(:,:,k) + dj;
    end
end

end

