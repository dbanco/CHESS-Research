function err = evalProbe(y,probe,N,D0,X,U)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
AD0 = reSampleNu(N,D0,probe(1),probe(2),U);
y_hat = ifft2(sum(bsxfun(@times,fft2(AD0),fft2(X)),3),'symmetric');
err = norm(y(:)-y_hat(:))/norm(y(:));
end

