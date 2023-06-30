function term = HornSchunckTermPaper(u,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
avgKernel = [1/12 1/6 1/12; 1/6 -1 1/6; 1/12 1/6 1/12];
uLap = convn(u,avgKernel,'same');
vLap = convn(v,avgKernel,'same');
term = sum((uLap).^2,'all') + sum((vLap).^2,'all');
end