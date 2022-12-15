function Xud = reSampleNuTrans(X,c1,c2,U)
%Dud = reSampleNu(N,D,c1,c2,U)
% N - Length of data
% D - Dictionary of dimensions [N1 x 1 x K]
% c1 - Integer factor 1
% c2 - Integer factor 2
% U - Number of scales (odd integer)

% N1 - Length of base dictonary entries
[N1,N2,KU] = size(X);
K = KU/U;

scales = ones(2,U);
V = ((U-1)/2);
% downscalings
for i = 1:V
    expon = (V-i+1);
    scales(:,i) = [c2^expon;c1^expon];
end
% upscalings
for i = (V+2):U
    expon = (i-V-1);
    scales(:,i) = [c1^expon;c2^expon];
end

Xud = zeros(N1,N2/4,K);
Xarray = reshape(X,[N1,N2,U,K]);
for k = 1:K
    for u = 1:U
        X_in = Xarray(:,:,u,k);
        Xud_k = uSampleTrans(dSampleTrans(X_in,scales(2,u)),scales(1,u));
        NN = min(size(Xud_k,2),N2/4);
        Xud(:,1:NN,k) = Xud(:,1:NN,k) + Xud_k(:,1:NN);
    end
end