function Xu = upSample(X,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2,K,T] = size(X);

if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end
Xu = zeros(N1,N2*2^(numScales-1),K*numScales,T);
% Construct stacked decimation matrices

for k = 1:K
    Xu(:,1:N2,1+numScales*(k-1),:) = X(:,:,k,:);
    for j = 2:numScales
        i = j+(k-1)*numScales;
        jj = 2^(j-1);
        Xu(:,2:2:N2*jj,i,:) = Xu(:,1:N2*(jj/2),i-1,:);
        Xu(:,1:N2*jj,i,:) = lwpass4(Xu(:,1:N2*jj,i,:),0);
    end
end

end

