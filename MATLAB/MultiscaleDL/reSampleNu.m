function Dout = reSampleNu(N,D,c1,c2,U)
%Dud = reSampleNu(N,D,c1,c2,U)
% N - Length of data
% D - Dictionary of dimensions [N1 x 1 x K]
% c1 - Integer factor 1
% c2 - Integer factor 2
% U - Number of scales (odd integer)

% N1 - Length of base dictonary entries
V = ((U-1)/2);
[N1,N2,K] = size(D);
Dud = zeros(N1,N,U,K);

for k = 1:K
    % downscalings
    for i = 1:V
        D2 = D(:,:,k);
        for ii = (V+1-i):-1:1
            D1 = uSample(D2,c2);
            D2 = dSample(D1,c1);
        end
        Nud = size(D2,2);
        Dud(:,1:Nud,i,k) = D2;
    end

    % neutral
    Dud(:,1:N2,V+1,k) = D(:,:,k);

    % upscalings
    for i = (V+2):U
        D2 = D(:,:,k);
        for ii = 1:(i-V-1)
            D1 = uSample(D2,c1);
            D2 = dSample(D1,c2);
        end
        Nud = size(D2,2);
        Dud(:,1:Nud,i,k) = D2;
    end
end

Dout = reshape(Dud,[N1,N,K*U]);




