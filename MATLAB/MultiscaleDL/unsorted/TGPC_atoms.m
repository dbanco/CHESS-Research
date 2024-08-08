%% Get learned dictionary atoms
% Gaussian peak functions with Poisson statistics

% T = 10;
% N = 101;
% [B,B_poiss,awmv_true] = generateExampleData(N,T);
% 
% figure(1)
% subplot(2,1,1)
% waterfall(B')
% subplot(2,1,2)
% waterfall(B_poiss')


topDir = 'C:\Users\dpqb1\Documents\Outputs2024\';
outDir = 'gaus_example_1_29_24_X0_D0_V00_sig_2';
folderPath = fullfile(topDir,outDir);
files = dir(fullfile(folderPath, '*.mat'));
matFileNames = {files.name};

load(fullfile(folderPath,matFileNames{4}))

D = outputs.D;
N = outputs.N;
M = outputs.M;
T = outputs.T;
y = outputs.y;
X = outputs.X;
K = outputs.K;
opt = outputs.opt;
scales = outputs.scales;
lambda = outputs.lambda;
lambda2 = outputs.lambda2;
Uvel = outputs.Uvel;
Vvel = outputs.Vvel;
center = (M+1)/2;

AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
Yhat = gather(Yhat);
err = sum((squeeze(y)-Yhat).^2,'all');

B(:,1) = squeeze(gather(D(1,:,:)));
B(:,2) = squeeze(gather(D(1,:,:)));
T = 2;
%% Define parameters

% Length of intensity data (theta coordinate)
P.N = size(B,1); 
N = P.N;

% Define dictionary of Gaussian basis functions
P.K = 20;   % Number of different basis functions 
P.sigmas = linspace(1/2,20,P.K); % Sigmas of basis functions

% ADMM parameters
params.adaptRho = 1; % binary flag for adaptive rho
params.mu = 2;       % tolerated factor between primal/dual residual
params.tau = 1.05;   % rho update factor
params.alpha = 1.8; % over-relaxation paramter
params.isNonnegative = 1; % flag to enforce nonnegativity
params.normData = 1; % flag to normalize b(t)

params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 100;

% Conjugate gradient parameters
params.conjGradIter = 100;
params.tolerance = 1e-8;
params.cgEpsilon = 1e-6;

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;

% Construct dictionary
A0ft = peakDictionaryFFT(P);
A0 = peakDictionary(P);

%% Setup and solve
% Without temporal coupling
params.rho1 = 1;                % ADMM parameter 1
params.lambda = 5e-2*ones(T,1); % Sparsity parameters
params.rho2 = 0;                % ADMM parameter 2
params.gamma = 0;               % Smoothness parameter
X_hat_indep = convADMM_LASSO_CG_TVphi_1D(A0ft,B,zeros(N,P.K,T),params);
B_hat_indep = Ax_ft_1D_Time(A0ft,X_hat_indep);

% With temporal coupling
% params.rho2 = 0.1;     % ADMM parameter 2
% params.gamma = 5e-3;% Smoothness parameter
% X_hat = convADMM_LASSO_CG_TVphi_1D(A0ft,B_poiss,zeros(N,P.K,T),params);
% B_hat = Ax_ft_1D_Time(A0ft,X_hat);

%% Get AWMV of atoms
awmv1 = computeAWMV_1D(X_hat_indep,P.sigmas);

%Compute AWMV at each scale
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],10,32, 1/10);
J = size(scales{1},2);

factors = scales{1}(1,:)./scales{1}(2,:);
awmvs = awmv1(1)*factors;


X_MCDL = gather(squeeze(X));
awmv_data = computeAWMV_1D(X_MCDL,awmvs);


%% Use atom AWMV to compute AWMV of full multiscale dictionary
position = 15+[linspace(20,50,20),30+20.*linspace(1,0.1,30).*cos((0:29)/3)];
true_std = (70-fliplr(position))/3;

% figure
% plot(true_std,'o-')

awmv_data = awmv_data-min(awmv_data);
awmv_data = awmv_data*25/max(awmv_data);
awmv_data = awmv_data+1;
hold on
plot(awmv_data,'-w','Linewidth',1)
% 
% fig3 = figure;
% plot(B(:,1))
% hold on
% plot(B_hat_indep(:,1))

%% Plot awmv recovery
% fig3 = figure;
% subplot(1,2,1)
% waterfall(B')
% title('truth')
% 
% subplot(1,2,2)
% waterfall(B_hat_indep')
% 
% 
% fig4 = figure;
% awmv1 = computeAWMV_1D(X_hat_indep,P.sigmas);
% plot(awmv1,'o-','Linewidth',2)
% ylabel('AWMV')
% xlabel('t')
% title('AWMV')
