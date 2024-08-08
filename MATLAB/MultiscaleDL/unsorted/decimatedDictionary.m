function [decimD,A] = decimatedDictionary(D,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

[N1,N2,K] = size(D);
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end

decimD = zeros(N1,N2,K*numScales);

% Construct stacked decimation matrices
A = eye(N2);
for k = 1:K
    for j = 1:numScales
        jj = 2^(j-1);
        i = j+(k-1)*numScales;
        B = zeros(N2);
        B(1:N2/jj,:) = A(jj:jj:N2,:);
        decimD(:,:,i) = B*squeeze(D(:,:,k))';
    end
end

end

