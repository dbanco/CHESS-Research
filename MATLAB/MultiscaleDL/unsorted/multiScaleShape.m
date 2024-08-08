function ydu = multiScaleShape(N,U)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f = zeros(1,N);
x = 1:N-16;
f(9:N-8) = exp(-x/(N))+0.2*sin(x/5);
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

yd = decimate(reshape(Pnrm(f),[1,N,1,1]),U);
ff = Pnrm(yd(:,1:N/2^(U-1),U));
ydu = upDwnSample(ff,U-1);
end

