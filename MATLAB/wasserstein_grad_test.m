%% Wasserstein gradient test
nn = 5;
h = repmat(1:nn,[nn,1]);
y = ones(nn,nn)/nn*nn;
% y = abs(randn(5,5)+0.5);
h = h/sum(h(:));
y = y/sum(y(:));


% Construct distance matrix
N = nn*nn;
THRESHOLD = nn*nn;
D = ones(N,N).*THRESHOLD;
for i = 1:nn
    for j = 1:nn
        for ii=max([1 i-THRESHOLD+1]):min([nn i+THRESHOLD-1])
            for jj = max([1 j-THRESHOLD+1]):min([nn j+THRESHOLD-1])
                ind1 = i + (j-1)*nn;
                ind2 = ii + (jj-1)*nn;
                D(ind1,ind2)= sqrt((i-ii)^2+(j-jj)^2); 
            end
        end
    end
end

deltaX = 1e-8;
lam = 10;
dFdX = zeros(size(h));

% [ gradW ] = WassersteinGrad( h(:), y(:), lam, D );
% [grad_a,grad_b] = WassersteinGrad2(h(:),y(:),lam,D);
% [gradFD,~] = WassersteinGradFD(h(:),y(:),lam,D);

% gradW = reshape(gradW,size(y));
% grad_a = reshape(grad_a,size(y))
% grad_b = reshape(grad_b,size(y))
% gradW
% gradFD*max(gradW(:))/max(gradFD(:))
lam = 0.01;
beta = 1;
[Wd_ipot,~,~,T1,a_t] = InexactProxOT(h,y,D,beta,1);
[Wd,~,~,T2] = sinkhornKnoppTransport(h,y,lam,D);
figure(11)
imagesc(T1)
title('IPOT')

figure(22)
imagesc(T2)
title('Sinkhorn')

[ gradW ] = WassersteinGrad( h(:), y(:), lam, D );
alpha_t = {};
gradIPOT = 0;
% for i = 1:numel(a_t)
%   gradIPOT = gradIPOT + (log(a_t{i})-log(sum(a_t{i})./numel(h)))/beta;
% end
% gradIPOT = gradIPOT/numel(a_t);
gradIPOT = gradIPOT + (log(a_t{end})-log(sum(a_t{end})./numel(y)))/beta;
figure(33)
imagesc(gradIPOT)
title('gradIPOT')

figure(44)
imagesc(gradW)
title('gradW')

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