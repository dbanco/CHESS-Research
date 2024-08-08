function Xd = decimate(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,K,T] = size(X);

if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xd = zeros(N1,N2,K*numScales,T);
% Construct stacked decimation matrices
for k = 1:K
    Xd(:,1:N2,1+numScales*(k-1),:) = X(:,:,k,:);
    for j = 2:numScales
        i = j+(k-1)*numScales;
        jj = 2^(j-1);
        xFilt = lwpass4(Xd(:,1:N2/jj*2,i-1,:),0);
        Xd(:,1:N2/jj,i,:) = xFilt(:,2:2:N2/jj*2,:);
    end
end

end

