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
Ft = -timeDiff(data);

u = zeros(size(data));
v = zeros(size(data));

avgKernel = [0 1 0; 1 0 1; 0 1 0]/4;
% avgKernel = [1 2 1; 2 -12 2; 1 2 1]/12;

for i = 1:maxIters
    % Update neighborhood averages

%     uPad = padarray(u,[1,1],'symmetric');
%     vPad = padarray(v,[1,1],'symmetric');
%     uAvg = convn(uPad,avgKernel,'valid');
%     vAvg = convn(vPad,avgKernel,'valid');

    uAvg = convn(u,avgKernel,'same');
    vAvg = convn(v,avgKernel,'same');
    factor = (Fx.*uAvg + Fy.*vAvg + Ft)./(smoothness + Fx.^2 + Fy.^2);
    factor(isnan(factor)) = 0;
    u = uAvg - Fx.*factor;
    v = vAvg - Fy.*factor;
end

end