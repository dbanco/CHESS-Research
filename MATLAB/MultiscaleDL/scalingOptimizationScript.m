%% Multiscale 1D dictionary learning toy problem
[y,yd,X_true,N,T,K,U,c1,c2] = rand_multirate_problem;
y = reshape(y,[1,N,1,T]); 
% plotDictUsage(yd,K,X_true)
% animateDataSeq(y)
V = (U-1)/2;
N0 = N/2^V;
D0 = zeros(1,N0,K);
D0(1,:,1) = yd(1,1:N0,V+1);% + 100*rand(1,N0,1)/100;
D0(1,:,2) = yd(1,1:N0,V+1+U);% + 100*rand(1,N0,1)/100;

% Compute y_hat
AD0 = reSampleNu(N,D0,c1,c2,U);
% plotDictUsage(AD0,K,X_true)
y_hat = ifft2(sum(bsxfun(@times,fft2(AD0),fft2(X_true)),3),'symmetric');

% Confirm reconstruction
norm(y(:)-y_hat(:))/norm(y(:))
% figure
% subplot(2,1,1)
% imagesc(squeeze(y))
% subplot(2,1,2)
% imagesc(squeeze(y_hat))
% 
% figure
% hold on
% plot(y(1,:,1,29),'Linewidth',2)
% plot(y_hat(1,:,1,29))

%%  Input: y,D0,N,U   Find c1,c2
denLim = 22;
% [output,err] = fareyScaleSearch(y,N,D0,X_true,U,denLim)
% [output] = goldenScaleSearch(y,N,D0,X_true,U,denLim)
%% Test how often we recover true solution
Trials = 1000;
trialErr = zeros(Trials,1);
for t = 1:Trials
    [y,yd,X_true,N,T,K,U,c1,c2] = rand_multirate_problem;
    denLim = 50;
    [output] = fareyScaleSearch(y,N,D0,X_true,U,denLim);
    trialErr(t) = norm(output - [c1;c2]);
end
plot(trialErr)
% close all
% Set up cbpdndl parameters
% lambda = 1e-1;
% K = 2;
% N0 = N/2^U;
% opt.DictFilterSizes = [ones(1,K);
%                        N0*ones(1,K)];
% opt.numScales = U;
% 
% 
% 
% plotDictUsage(Dr,K,X_true)