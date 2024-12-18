function [i1,j1,i2,j2] = getMin(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
data(data == 0) = nan;
d1 = data(:,:,1);
d2 = data(:,:,2);
[N1,N2] = size(d1);
ind1 = find(data(:,:,1) == min(d1(:)));
ind2 = find(data(:,:,2) == min(d2(:)));
[i1,j1] = ind2sub([N1,N2],ind1);
[i2,j2] = ind2sub([N1,N2],ind2);
end