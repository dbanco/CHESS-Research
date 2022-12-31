function Xout = reSampleNuTrans2(M,X,c1,c2,U,NormVals)
%Dud = reSampleNu(N,D,c1,c2,U)
% N - Length of data
% D - Dictionary of dimensions [N1 x 1 x K]
% c1 - Integer factor 1
% c2 - Integer factor 2
% U - Number of scales (odd integer)

% N1 - Length of base dictonary entries


V = ((U-1)/2);
[N1,N2,KU] = size(X);
K = KU/U;
% N = floor(N2*c1^V/c2^V);
Xud = zeros(N1,N2,K);
Xarray = reshape(X,[N1,N2,U,K]);

for k = 1:K
    % neutral
    Xud(:,1:N2,k) = Xarray(:,1:N2,V+1,k);
    
    % downscalings
    for i = 1:V
        Xin = Xarray(:,:,i,k);
        for ii = (V+1-i):-1:1
            Xin2 = dSampleTrans(M,Xin,c1);
            clear Xin
            Xin = uSampleTrans(M,Xin2,c2);
            clear Xin2
        end
        Nud = min(size(Xin,2),N2);
        Xud(:,1:Nud,k) = Xud(:,1:Nud,k) + Xin(:,1:Nud)/NormVals(i+U*(k-1));
    end

    % upscalings
    for i = (V+2):U
        Xin = Xarray(:,:,i,k);
        for ii = 1:(i-V-1)
            Xin2 = dSampleTrans(M,Xin,c2);
            clear Xin
            Xin = uSampleTrans(M,Xin2,c1);
            clear Xin2
        end
        Nud = min(size(Xin,2),N2);
        Xud(:,1:Nud,k) = Xud(:,1:Nud,k) + Xin(:,1:Nud)/NormVals(i+U*(k-1));
    end
end
Xout = Xud(N1,1:M,:);




