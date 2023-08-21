function Xu = uSampleTransCenter(M,X,m,center)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2] = size(X);
if m == 1
    Xu = X;
else
    Xu = zeros(N1,N2);
    d1 = lowpassM(X,m,1);
    inds = [fliplr(center-m:-m:1),(center:m:N2)];
    Xu = d1(inds);
end


end

