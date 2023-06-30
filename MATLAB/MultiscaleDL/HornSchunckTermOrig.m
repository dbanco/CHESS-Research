function term = HornSchunckTermOrig(u,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
avgKernel = [0 1 0; 1 0 1; 0 1 0]/4;
uAvg = convn(u,avgKernel,'same');
vAvg = convn(v,avgKernel,'same');
term = sum((uAvg-u).^2,'all') + sum((vAvg-v).^2,'all');
end