function dd = dSample(d,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

if m == 1
    dd = d;
else
    d1 = lowpassM(d,m,0);
    dd = d1(m:m:end);
end

end

