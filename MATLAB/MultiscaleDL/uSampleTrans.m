function Xu = uSampleTrans(M,X,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2] = size(X);
Xu = zeros(N1,N2);
d1 = lowpassM(X,m,1);
Xu = d1(m:m:N2);


end

