function Xd = dSampleTrans(M,X,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2] = size(X);
Xd = zeros(N1,m*N2);

Xd(:,m*(1:N2)) = X(:,1:N2);
Xd = lowpassM(Xd,m,1);

end

