function waveletScale = computeWaveletScale(Xout)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
N = size(Xout,2);
waveletScale = sum(sum(abs(Xout),1).*(N:-1:1))/sum(abs(Xout),'all');

end

