function f = multiRateShape(M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f = zeros(1,M);
x = 1:M;
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

f = Pnrm(5*M/4 - x + M/4*sin(x));

end

