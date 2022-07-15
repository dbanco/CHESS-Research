function ydu = multiRateSinc(N,scales)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f = zeros(1,N);
x = 1:N-16;
f(9:N-8) = sinc(8*x/N - pi) + 0.25;
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
fd = reSample(N,Pnrm(f(1:N/2)),[1;2]);
ydu = reSample(N,Pnrm(fd),scales);
end

