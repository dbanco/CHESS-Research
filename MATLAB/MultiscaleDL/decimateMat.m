function A = decimateMat(N,numScales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

if ceil(log2(N)) ~= floor(log2(N))
    error('Data size is not power of 2')
end

% Construct stacked decimation matrices
A = zeros(N,N,numScales);
A = eye(N);
for j = 1:numScales
    jj = 2^(j-1);
    A(1:N/jj,:,j) = A(jj:jj:N,:,1);
end


end

