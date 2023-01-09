function Xout = reSampleNuTrans2d(M,X,c1,c2,U,NormVals)
%Dud = reSampleNu(N,D,c1,c2,U)
% N - Length of data
% D - Dictionary of dimensions [N1 x 1 x K]
% c1 - Integer factor 1
% c2 - Integer factor 2
% U - Number of scales (odd integer)

% N1 - Length of base dictonary entries


V = ((U-1)/2);
[N1,N2,KU] = size(X);
N3 = max(N1,N2);
if N1>N2
    X = reshape(X,[N2,N1,KU]);
end


K = KU/U;
% N = floor(N2*c1^V/c2^V);
Xud = zeros(1,N3,K);
Xarray = reshape(X,[1,N3,U,K]);

for k = 1:K
    % neutral
    Xud(:,:,k) = Xarray(:,:,V+1,k);
    
    % downscalings
    for i = 1:V
        Xin = Xarray(:,:,i,k);
        for ii = (V+1-i):-1:1
            Xin2 = dSampleTrans(M,Xin,c1);
            clear Xin
            Xin = uSampleTrans(M,Xin2,c2);
            clear Xin2
        end
        Nud = min(size(Xin,2),N3);
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
        Nud = min(size(Xin,2),N3);
        Xud(:,1:Nud,k) = Xud(:,1:Nud,k) + Xin(:,1:Nud)/NormVals(i+U*(k-1));
    end
end
Xout = Xud(1,1:M,:);

if N1>N2
    Xout = reshape(Xout,[M,1,K]);
end



