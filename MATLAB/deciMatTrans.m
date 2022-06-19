function Atv = deciMatTrans(A,v)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

[N1,N2,KU,T] = size(v);
numScales = size(A,3);
K = KU/numScales;
if ceil(log2(N2)) ~= floor(log2(N2))
    error('Data size is not power of 2')
end

Atv = zeros(N1,N2,K,T);

% Construct stacked decimation matrices
for k = 1:K
    for j = 1:numScales
        i = j+(k-1)*numScales;
        At = squeeze(A(:,:,j))';
        c = reshape(v(:,:,i,:),[N2,T]);
        Atv(:,:,k,:) = Atv(:,:,k,:) + reshape(At*c,[N1,N2,1,T]);    
    end
end

end