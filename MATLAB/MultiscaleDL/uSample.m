function du = uSample(d,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here
[N1,N2] = size(d);

du = zeros(N1,N2*m);

du(:,m*(1:N2)) = d;
du = lowpassM(du,m,1);


end

