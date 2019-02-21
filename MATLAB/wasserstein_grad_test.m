%% Wasserstein gradient test
h = repmat(1:5,[5,1]);
y = ones(5,5)/25;
% y = abs(randn(5,5)+0.5);
h = h/sum(h(:));
y = y/sum(y(:));


% Construct distance matrix
N = 25;
THRESHOLD = 25;
D = ones(N,N).*THRESHOLD;
for i = 1:5
    for j = 1:5
        for ii=max([1 i-THRESHOLD+1]):min([5 i+THRESHOLD-1])
            for jj = max([1 j-THRESHOLD+1]):min([5 j+THRESHOLD-1])
                ind1 = i + (j-1)*5;
                ind2 = ii + (jj-1)*5;
                D(ind1,ind2)= sqrt((i-ii)^2+(j-jj)^2); 
            end
        end
    end
end

deltaX = 1e-8;
lam = 10;
dFdX = zeros(size(h));

[ gradW ] = WassersteinGrad( h(:), y(:), lam, D );
[grad_a,grad_b] = WassersteinGrad2(h(:),y(:),lam,D);
[gradFD,~] = WassersteinGradFD(h(:),y(:),lam,D);

gradW = reshape(gradW,size(y));
grad_a = reshape(grad_a,size(y))
grad_b = reshape(grad_b,size(y))
gradW
gradFD*max(gradW(:))/max(gradFD(:))
%%
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