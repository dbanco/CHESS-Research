function Xout = reSampleTransCustomArrayCenter2(M,X,scales,centers,NormVals)
%Dud = reSampleNu(N,D,c1,c2,U)
% M - Length of base dictionary entries
% X - Dictionary of dimensions [N1 x N2 x KJ]
% scales - Array of scaling fators [2 x U]
% NormVals (optional) - Array of normalization factors [U x 1]

[N1,N2,Utotal] = size(X);
K = numel(scales);
Uarray = zeros(K,1);
for i = 1:K
    Uarray(i) = size(scales{i},2);
end

if nargin < 5
    warning('Normalization factors not provided')
    NormVals = ones(Utotal,1);
end

Xud = zeros(N1,N2,K);

ind1 = 1;
for k = 1:K
    U = Uarray(k);
    ind2 = ind1 + U - 1;
    scales_k = scales{k};
    Xarray_k = X(:,:,ind1:ind2);

    for i = 1:U
        % Get numerator and denominator
        c1 = scales_k(1,i);
        c2 = scales_k(2,i);

        % Apply rescalings
        if c1 == c2
            Xud(:,:,k) = Xud(:,:,k) + Xarray_k(:,:,i);
        else

            Xin = Xarray_k(:,:,i);
            [Xin2,center_up] = uSampleCenter2(Xin,c2,c2*centers(i)); clear Xin
            [Xin,center_dwn] = dSampleCenter2(Xin2,c1,center_up); clear Xin2
            shift = center_dwn-centers(i);
            Xud(:,:,k) = Xud(:,:,k) + Xin(:,1+shift:N2+shift)/NormVals(i+ind1-1);
        end
    end
    ind1 = ind1 + Uarray(k);
end
Xout = real(Xud(:,1:M,:));




