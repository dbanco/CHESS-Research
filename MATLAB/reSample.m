function Dud = reSample(N,D,scales)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,~,K] = size(D);
U = size(scales,2);
Dud = zeros(N1,N,U,K);

for k = 1:K
    for u = 1:U
        D1 = uSample(D(:,:,k),scales(1,u));
        D2 = dSample(D1,scales(2,u));
        Nud = min(N,numel(D2));
        Dud(:,1:Nud,u,k) = D2(1:Nud);
    end
end

Dud = reshape(Dud,[N1,N,K*U]);




