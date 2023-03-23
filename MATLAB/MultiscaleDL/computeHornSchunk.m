function [u,v,Fy,Fx,Ft] = computeHornSchunk(data,smoothness,maxIters)
%computeHornSchunk Summary of this function goes here
%   Detailed explanation goes here

if nargin <2
    smoothness = 1;
end
if nargin <3
    maxIters = 10;
end

Fy = sobel(data,1);
Fx = sobel(data,2);
Ft = timeDiff(data);

u = zeros(size(data));
v = zeros(size(data));
avgKernel = [0 1 0; 1 0 1; 0 1 0];

for i = 1:maxIters
    % Update neighborhood averages
    uAvg = convn(u,avgKernel,'same');
    vAvg = convn(v,avgKernel,'same');

    u = uAvg - Fx.*(Fx.*uAvg + Fy.*vAvg + Ft)./(smoothness^2 + Fx.^2 + Fy.^2);
    v = vAvg - Fy.*(Fx.*uAvg + Fy.*vAvg + Ft)./(smoothness^2 + Fx.^2 + Fy.^2);
end

end