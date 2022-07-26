function [P,N,K,M,T] = definePoissonP(P)
% File Parameters
P.baseFileName = 'indep_fit_%i_%i.mat';
num_ims = 50;
N = 101;
K = 20;
M = 50;
T = num_ims;
zPad = 0;
zMask = [];
P.dataScale = 1;
P.lambda_values = logspace(-3,1,M);
P.num_theta = N;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = linspace(0.5,25,P.num_var_t).^2;

% algorithm parameters

P.params.rho1 = 1;
% P.params.lambda1 = 0.0001;
P.params.tau = 1.05;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 50;
P.params.tolerance = 1e-8;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

P.params.conjGradIter = 100;
P.params.tolerance = 1e-8;
P.params.cgEpsilon = 1e-3;
end

