function AD = deciMat(A,D)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

[N1,N2,K] = size(D);
numScales = size(A,3);
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end

AD = zeros(N1,N2,K*numScales);

% Construct stacked decimation matrices
for k = 1:K
    for j = 1:numScales
        i = j+(k-1)*numScales;
        AD(:,:,i) = A(:,:,j)*squeeze(D(:,:,k))';
    end
end

end