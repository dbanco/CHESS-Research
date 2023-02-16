function Xout = reSampleTransCustom(M,X,scales,NormVals)
%Dud = reSampleNu(N,D,c1,c2,U)
% M - Length of base dictionary entries
% D - Dictionary of dimensions [N1 x M x K]
% scales - Array of scaling fators [2 x U]
% NormVals (optional) - Array of normalization factors [U x 1]

U = size(scales,2);
[N1,N2,KU] = size(X);
K = KU/U;

if nargin < 4
    warning('Normalization factors not provided')
    NormVals = ones(K*U,1);
end

Xud = zeros(N1,N2,K);
Xarray = reshape(X,[N1,N2,U,K]);

for k = 1:K
    % downscalings
    for i = 1:U
        c1 = scales(1,i);
        c2 = scales(2,i);
        if c1 == c2
            Xud(:,:,k) = Xud(:,:,k) + Xarray(:,:,i,k);
        else
            Xin = Xarray(:,:,i,k);
            Xin2 = dSampleTrans(M,Xin,c2); clear Xin
            Xin = uSampleTrans(M,Xin2,c1); clear Xin2
            Nud = min(size(Xin,2),N2);
            Xud(:,1:Nud,k) = Xud(:,1:Nud,k) + Xin(:,1:Nud)/NormVals(i+U*(k-1));
        end
    end
end
Xout = Xud(N1,1:M,:);




