function y = filterbankSyn(Xout,bank)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[h1,h2] = filters(bank);
N = size(h1,1);


M = size(Xout,1);
xLoLo = Xout(:,1);
xLoHi = Xout(:,2);
xHiLo = Xout(:,3);
xHiHi = Xout(:,4);

xLo = conv(upsample(xLoLo,2),h1(:,2)) + ...
      conv(upsample(xLoHi,2),h2(:,2));
xHi = conv(upsample(xHiLo,2),h1(:,2)) + ...
      conv(upsample(xHiHi,2),h2(:,2));

y = conv(upsample(xLo,2),h1(:,2)) + ...
    conv(upsample(xHi,2),h2(:,2));
y = y(1+3*(N-1):end-3*(N-1));
end
