function [Dout,NormVals] = reSampleCustom(N,D,scales,NormVals)
%Dud = reSampleNu(N,D,c1,c2,U)
% N - Length of data
% D - Dictionary of dimensions [N1 x 1 x K]
% scales - Array of scaling fators [2 x U]
% NormVals (optional) - Array of normalization factors [U x 1]

[N1,N2,K] = size(D);
U = size(scales,2);
Dout = zeros(N1,N,U,K);


for k = 1:K
    % downscalings
    for i = 1:U
        c1 = scales(1,i);
        c2 = scales(2,i);
        if c1 == c2
            Dout(1:N1,1:N2,i,k) = D(1:N1,1:N2,k);
        else
            Dk = D(:,:,k);
            Dku = uSample(Dk,c1);
            Dk = dSample(Dku,c2);
            [Nud1,Nud2,~] = size(Dk);
            Dout(1:Nud1,1:Nud2,i,k) = Dk;
        end
    end
end

Dout = reshape(Dout,[N1,N,K*U]);

if nargin ==4
    for i = 1:K*U
        Dout(:,:,i) = Dout(:,:,i)/NormVals(i);
    end
else
    NormVals = zeros(K*U,1);
    for i = 1:K*U
        NormVals(i) = norm(Dout(:,:,i));
        Dout(:,:,i) = Dout(:,:,i)/NormVals(i);
    end
end




