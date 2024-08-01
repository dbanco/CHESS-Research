function ydu = multScaleShape(N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f(9:N-8) = sinc(8*x/N - pi) + 0.25;
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

yd = decimate(reshape(Pnrm(f),[1,N,1,1]),5);
ff = Pnrm(yd(:,1:N/4,3));
ydu = upDwnSample(ff,2);
end

