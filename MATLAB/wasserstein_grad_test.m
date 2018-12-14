%% Wasserstein gradient test
h = repmat(1:5,[5,1]);
y = ones(5,5)/25;
% y = abs(randn(5,5)+0.5);
h = h/sum(h(:));
y = y/sum(y(:));

lam = 1;
% Construct distance matrix
THRESHOLD = 32;
N = 25;
D = ones(N,N).*THRESHOLD;
for i=1:N
    for j=max([1 i-THRESHOLD+1]):min([N i+THRESHOLD-1])
        D(i,j)= abs(i-j); 
    end
end
% D = D/max(D(:)+0.1);

[ gradW ] = WassersteinGrad( y(:), h(:), lam, D );
gradW = gradW/sum(abs(gradW));


figure(1)
subplot(2,2,1)
imagesc(h)
colorbar()
title('P_1')
subplot(2,2,2)
imagesc(y)
colorbar()
title('P_2')
subplot(2,2,3)
imagesc(reshape(gradW,[5,5]))
colorbar()
title(sprintf('grad W for lambda=%0.000f',lam))
subplot(2,2,4)
imagesc(h+reshape(gradW,[5,5]))
colorbar()
title('P_1 - gradW')