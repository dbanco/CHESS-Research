function ydu = multiRateShapeNu(N,c1,c2,U)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f = zeros(1,N);
x = 1:N-16;
f(9:N-8) = exp(-x/(N))+0.2*sin(x/5);
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
fd = reSample(N,Pnrm(f(1:N/2)),[1;2]);
ydu = reSampleNu(N,Pnrm(fd(1:N/4)),c1,c2,U);
end

