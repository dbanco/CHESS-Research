function [u,v,Fx,Fy,Ft] = computeHornSchunkPaper(data,smoothness,maxIters)
%computeHornSchunk Summary of this function goes here
%   Detailed explanation goes here

if nargin <2
    smoothness = 1;
end
if nargin <3
    maxIters = 10;
end

dataPad = padarray(data,[1 1 1],0,'pre');
Fx = diffxHS(dataPad);
Fy = diffyHS(dataPad);
Ft = difftHS(dataPad);

u = zeros(size(data));
v = zeros(size(data));

avgKernel = [1 2 1; 2 0 2; 1 2 1]/12;

for i = 1:maxIters
    % Update neighborhood averages
    uAvg = convn(u,avgKernel,'same');
    vAvg = convn(v,avgKernel,'same');
    
    factor = (Fx.*uAvg + Fy.*vAvg + Ft)./(smoothness + Fx.^2 + Fy.^2);
    factor(isnan(factor)) = 0;
    u = uAvg - Fx.*factor;
    v = vAvg - Fy.*factor;

    
%     [obj1, obj2, sys] = HSobjectivePaper(Fx,Fy,Ft,u,v,1,smoothness);
%     ofObj = obj1+obj2*smoothness
%     hsObj = norm(sys(:))
end

end