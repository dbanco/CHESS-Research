function out = AfuncTrans(in,U,xf,N,T)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

in = reshape(in(:),[1,N,1,T]);
out2 = decimateTrans(ifft2(sum(bsxfun(@times, xf, fft2(in)), 4),'symmetric'),U);
out = out2(:);


