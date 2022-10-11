%% Find scale parameter golden section search

%% Multiscale 1D dictionary learning toy problem
[y,yd,X_true,N,T,K,c1,c2] = rand_multirate_problem;
y = reshape(y,[1,N,1,T]); 
plotDictUsage(yd,K,X_true)
% animateDataSeq(y)

%%
close all
% Set up cbpdndl parameters
lambda = 1e-1;
K = 2;
N0 = N/2^U;
opt.DictFilterSizes = [ones(1,K);
                       N0*ones(1,K)];
opt.numScales = U;

D0 = zeros(1,N0,K);
D0(1,:,1) = yd(1,1:N0,U+1);% + 100*rand(1,N0,1)/100;
D0(1,:,2) = yd(1,1:N0,U+1+2*U);% + 100*rand(1,N0,1)/100;

sampFactors = [1 1 1 2;...
               4 2 1 1];
Dr = reSample(N,D0,sampFactors);
plotDictUsage(Dr,K,X_true)