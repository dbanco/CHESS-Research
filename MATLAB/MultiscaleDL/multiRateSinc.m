function f = multiRateSinc(M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f = zeros(1,M);
x = 0.5:M-0.5;
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

f = Pnrm(sinc(2*pi*x/M - pi) + 0.25);


end

