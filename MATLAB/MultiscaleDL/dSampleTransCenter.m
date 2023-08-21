function Xd = dSampleTransCenter(M,X,m,center)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2] = size(X);
Xd = zeros(N1,m*N2);
inds = [fliplr(center-m:-m:1),center:m:m*N2];
Xd(:,inds) = X(:,1:N2);
Xd = lowpassM(Xd,m,1);

end

