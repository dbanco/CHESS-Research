function checkSolution(y,X,D,U,f)
AD0 = decimate(D,U);
AD0f = fft2(AD0);
Yhat0 = squeeze(ifft2(sum(bsxfun(@times,AD0f,fft2(X)),3),'symmetric'));

if nargin <5
   f = figure;
else
   f = figure(f);
end
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat0)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))
p = 1+(f.Number-1)*400;
f.Position = [p 1 400 500];
end

