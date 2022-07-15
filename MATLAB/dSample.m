function dd = dSample(d,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

d1 = lowpassM(d,m,1);
dd = d1(1:m:end);

end

