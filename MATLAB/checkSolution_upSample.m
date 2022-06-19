function checkSolution_upSample(y,X,D,U)
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
AD0 = upSample(D,U);
AD0f = fft2(AD0);
Yhat0 = squeeze(ifft2(sum(bsxfun(@times,AD0f,fft2(X)),3),'symmetric'));

figure(10);
subplot(1,2,1)
imagesc(squeeze(y))
subplot(1,2,2)
imagesc(Yhat0)
title(sprintf('Rel Error: %0.3f',norm(squeeze(y)-Yhat0,'fro')/norm(y(:),'fro')))

end

