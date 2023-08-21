function [Dout,NormVals] = reSampleCustomArray(N,D,scales,NormVals)
%Dud = reSampleCustomArray(N,D,scales,NormVals)
% N - Length of data
% D - Dictionary of dimensions [N1 x 1 x K]
% scales_k - Cell array of arrays of scaling fators [2 x Ui]
% NormVals (optional) - Array of normalization factors [U x 1]
if numel(size(D)) == 2
    N1 = 1;
    [N2,K] = size(D);
    D = reshape(D,[N1,N2,K]);
else
    [N1,N2,K] = size(D);
end

Uarray = zeros(K,1);
for i = 1:K
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);

Dout = zeros(N1,N,Utotal);

ind1 = 1;
for k = 1:K
    U = Uarray(k);
    ind2 = ind1 + U - 1;
    Dout_k = zeros(N1,N,U);
    scales_k = scales{k};

    for i = 1:U
        % Get numerator and denominator
        c1 = scales_k(1,i);
        c2 = scales_k(2,i);

        % Apply rescalings
        if c1 == c2
            Dout_k(1:N1,1:N2,i) = D(1:N1,1:N2,k);
        else
            Dk = D(:,:,k);
            Dku = uSample(Dk,c1);
            Dk = dSample(Dku,c2);
            [Nud1,Nud2,~] = size(Dk);
            Dout_k(1:Nud1,1:min(Nud2,N),i) = Dk(1:Nud1,1:min(Nud2,N));
        end
    end
    Dout(:,:,ind1:ind2) = Dout_k;
    ind1 = ind1 + U;
end

if nargin ==4
    for i = 1:Utotal
        Dout(:,:,i) = Dout(:,:,i)/NormVals(i);
    end
else
    NormVals = zeros(Utotal,1);
    for i = 1:Utotal
        NormVals(i) = norm(Dout(:,:,i));
        Dout(:,:,i) = Dout(:,:,i)/NormVals(i);
    end
end




