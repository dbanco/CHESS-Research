% Dictionary Constrained Nonnegative CP
% Load MMPAD data
data_dir = 'D:\MMPAD_data_nr1\ring1_zero';
t = 1:4:200;
T = numel(t);
for i = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t(i)),'.mat']))
    if i == 1
        [N,M] = size(polar_image(:,134:end));
        Y = zeros(N,M,T);
    end
    Y(:,:,i) = polar_image(:,134:end)/norm(polar_image(:,134:end),'fro');
end

% Parameters
param.R = 20;
param.maxIters = 1000;
param.stopCrit = 1e-8;
param.verbose = 1;

param.lambda = 0.1;
param.phi1 = 1;
param.phi2 = 1;
param.kappa = 0;

param.ncpIters = 1;
param.ncpTol = 5e-4;
param.ncpVerbose = 0;
param.Uinit = 'rand';

% Set up cbpdndl parameters
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 1;
opt.rho = 50*param.lambda + 0.5;
opt.sigma = size(Y,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.RelStopTol = 0.005;

opt1 = opt; K1=9;
opt1.DictFilterSizes = [ones(1,K1);
                       8,8,8,16,16,16,32,32,32];
opt2 = opt; K2=12;
opt2.DictFilterSizes = [ones(1,K2);
                        16,16,16,32,32,32,64,64,64,128,128,128];
param.opt1 = opt1;
param.opt2 = opt2;


%exp13 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1  0.1 0.1 0.1 0.1 0.1 0.1 ];
loop1.kappa =  [0.15 0.12 0.10 0.08 0.06 0.04 0.02 0.01 0.05 0.02 0.01 0];
loop1.phi1 =   [0 0 0 0 0 0 0 0 0 0 0 0];
loop1.phi2 =   [0 0 0 0 0 0 0 0 0 0 0 0];

saveOuts = 1;
outDir = 'exp13';
topDir = ['C:\Users\dpqb1\Desktop\cdc_ncp_outputs\',outDir];
    mkdir(topDir)
for ii = 9:numel(loop1.lambda)
    ending = ['_out_',num2str(ii)];
    
    param.lambda = loop1.lambda(ii);
    param.kappa = loop1.kappa(ii);
    param.phi1 = loop1.phi1(ii);
    param.phi2 = loop1.phi2(ii);
    % Gaussians
%     % Initial 2theta dictionary
%     D1 = initGaussianDictionary(param.opt1.DictFilterSizes);  
%     % Initial eta dictionary
%     D2 = initGaussianDictionary(param.opt2.DictFilterSizes);

    % Ones
%     % Initial 2theta dictionary
%     Nd1 = max(opt1.DictFilterSizes(:));
%     D1 = ones(1,Nd1,K1);
%     % Initial eta dictionary
%     Nd2 = max(opt2.DictFilterSizes(:));
%     D2 = ones(1,Nd2,K2);

    % Random Uniform
    % Initial 2theta dictionary
    Nd1 = max(opt1.DictFilterSizes(:));
    D1 = rand(1,Nd1,K1);
    % Initial eta dictionary
    Nd2 = max(opt2.DictFilterSizes(:));
    D2 = rand(1,Nd2,K2);
    
    %% Solve cdc_ncp
    [U, D1, D2, X1, X2, Vals] = cdc_ncp(Y,D1,D2,param);
%     load(fullfile(topDir,['variables',ending,'.mat']))
    if saveOuts
        saveOutputs(topDir,ending,Y,U,D1,D2,X1,X2,Vals,param)
    end
end

% History

%exp1 - varying phi
% loop1.phi1 = [0, logspace(-3,log10(5),9)];
% loop1.phi2 = [0, logspace(-3,log10(5),9)];
% loop1.kappa = 0;

%exp2 - 500 outer, 5 inner (this was just bad)
% loop1.phi1 = [0,1.7242];
% loop1.phi2 = loop1.phi1;
% loop1.kappa = 0;

%exp3 - 100 outer, 500 inner w/ tol
% loop1.phi1 = [0, 1.7242];
% loop1.phi2 = loop1.phi1;
% loop1.kappa = 0;

%exp5 - 100 outer, 500 inner w/ tol
% loop1.lambda = logspace(-2,-0.5,10);
% param.phi1 = 1.7242;
% param.phi2 = 1.7242;
% param.kappa = 0;

%exp6 - 100 outer, 500 inner w/ tol
% loop1.phi1 = [0, 1.7242];
% loop1.phi2 = loop1.phi1;
% loop1.kappa = 0;

%exp7
% loop1.lambda = ones(1,3)*0.1;
% loop1.kappa = [0, 0.5, 0.25];
% loop1.phi1 = ones(1,3)*1.7242;
% loop1.phi2 = ones(1,3)*1.7242;

%exp8 200,200,200
% loop1.lambda = [0 0.1 0.1];
% loop1.kappa = zeros(1,3);
% loop1.phi1 = [0 0.5 1];
% loop1.phi2 = [0 0.5 1];

%exp8b same but with ones dict init
% loop1.lambda = [0 0.1 0.1];
% loop1.kappa = zeros(1,3);
% loop1.phi1 = [0 0.5 1];
% loop1.phi2 = [0 0.5 1];

%exp8c same but with rand uniform dict init
% loop1.lambda = [0 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa = zeros(1,6);
% loop1.phi1 = [0 0.5 1 1.2 1.4 1.6];
% loop1.phi2 = [0 0.5 1 1.2 1.4 1.6];


%exp9 rand uniform dict init, 1max inner iters
% loop1.lambda = [0 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa = zeros(1,6);
% loop1.phi1 = [0 0.5 1 1.2 1.4 1.6];
% loop1.phi2 = [0 0.5 1 1.2 1.4 1.6];

%exp10 rand uniform dict init, 5max inner iters
% loop1.lambda = [0 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa = zeros(1,6);
% loop1.phi1 = [0 0.5 1 1.2 1.4 1.6];
% loop1.phi2 = [0 0.5 1 1.2 1.4 1.6];

%exp11 rand uniform dict init, 1max inner iters
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1];
% loop1.kappa = [0.01 0.05 0.1 0.2 0.5];
% loop1.phi1 = zeros(1,5);
% loop1.phi2 = zeros(1,5);

%exp12 rand uniform dict init, 1max inner iters,param.Uinit = 'ones';
% loop1.lambda = [0.1  0.1];
% loop1.kappa =  [0.15 0.14];
% loop1.phi1 =   [1.4  1.4];
% loop1.phi2 =   [1.4  1.4];

%exp12b rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% loop1.lambda = [0.1  0.1 0.1  0.1 0.1  0.1];
% loop1.kappa =  [0.15 0.14 0.13 0.12 0.11 0.1];
% loop1.phi1 =   [1.4  1.4 1.4  1.4 1.4  1.4];
% loop1.phi2 =   [1.4  1.4 1.4  1.4 1.4  1.4];
