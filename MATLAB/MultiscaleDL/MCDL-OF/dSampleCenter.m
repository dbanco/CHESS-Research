function [dd,center] = dSampleCenter(d,m,center)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

if m == 1
    dd = d;
else
    N = numel(d);
    inds = [fliplr(center-m:-m:1),center:m:N];
    center = find(inds==center);
    d1 = lowpassM(d,m,0);
    dd = d1(inds);
end

end

