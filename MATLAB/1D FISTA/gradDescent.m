function [xkp1] = GradDescent( A,b,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L = params.L;
maxGradIters = params.maxGradIters;
gradTolerance = params.gradTolerance;


c =A'*R;
xk = zeros(N,M);
xkp1 = zeros(N,M);
c = A'*b;
AA = A'*A;
for j = 1:50
%     gradient = A'*(A*xk-b);
    gradient = AA*xk - c
    xkp1 = xk - (1/L)*gradient;
end
end

