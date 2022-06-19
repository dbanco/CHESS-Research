% Dictionary Constrained Nonnegative CP
%Generate toy data
Y = genSpotSeq2noisenorm;

% Parameters
param.R = 20;
param.maxIters = 1000;
param.stopCrit = 1e-7;
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
opt.MaxMainIter = 0;
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

opt1 = opt; K1=6;
opt1.DictFilterSizes = [ones(1,K1);
                        ones(1,K1)*64];
opt2 = opt; K2=8;
opt2.DictFilterSizes = [ones(1,K2);
                        ones(1,K2)*256];
param.opt1 = opt1;
param.opt2 = opt2;

%exp7 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% opt.MaxMainIter = 0;
% outDir = 'exp7';
% param.R = 20;
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa =  [0   0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1   1.8 2   10  20  50  100 200 500 1000];
% loop1.phi1 =   [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 ];
% loop1.phi2 =   [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];

%exp8 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Do better init and many iterations
% Y = genSpotSeq2noisenorm;
% opt.MaxMainIter = 20; with init iters 20,50,50
outDir = 'exp9';
param.R = 20;
loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
loop1.kappa =  [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
loop1.phi1 =   [logspace(-2,3,20) 0];
loop1.phi2 =   [logspace(-2,3,20) 0];


saveOuts = 1;
topDir = ['C:\Users\dpqb1\Desktop\cdc_ncp_toy1\',outDir];
mkdir(topDir)

for ii = 21:numel(loop1.lambda)
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

%exp1 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa =  [0   0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
% loop1.phi1 =   [0   0   0   0   0   0   0   0   0   0 ];
% loop1.phi2 =   [0   0   0   0   0   0   0   0   0   0 ];

%exp2 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa =  [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
% loop1.phi1 =   [0   0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1   1.5 2   2.5 3   5   10  20  40];
% loop1.phi2 =   [0   0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1   1.5 2   2.5 3   5   10  20  40];

%exp3 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]*5;
% loop1.kappa =  [0   0   0   0   0   0   0   0];
% loop1.phi1 =   [1.5 2   2.5 3   5   10  20  40];
% loop1.phi2 =   [1.5 2   2.5 3   5   10  20  40];

%exp4 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]*2   ;
% loop1.kappa =  [0   0   0   0   0   0   0   0   0];
% loop1.phi1 =   [1.5 2   2.5 3   5   10  20  40  100];
% loop1.phi2 =   [1.5 2   2.5 3   5   10  20  40  100];

%exp5_1 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% outDir = 'exp5_1';
% param.R = 3;
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 ]*2   ;
% loop1.kappa =  [0   0   0.2 0.2 0.3 0.3 ];
% loop1.phi1 =   [1.5 3   1.5 3   1.5 3];
% loop1.phi2 =   [1.5 3   1.5 3   1.5 3];

%exp5_2 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% outDir = 'exp5_2';
% param.R = 20;
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 ]*2   ;
% loop1.kappa =  [0   0   0.2 0.2 0.3 0.3 ];
% loop1.phi1 =   [1.5 3   1.5 3   1.5 3];
% loop1.phi2 =   [1.5 3   1.5 3   1.5 3];

%exp5_3 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% outDir = 'exp5_3';
% param.R = 20;
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 ]*2   ;
% loop1.kappa =  [0   0   0.2 0.2 0.3 0.3 ];
% loop1.phi1 =   [1.5 3   1.5 3   1.5 3];
% loop1.phi2 =   [1.5 3   1.5 3   1.5 3];

%exp5_4 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% outDir = 'exp5_4';
% param.R = 20;
% loop1.lambda = [0.1 0.1];
% loop1.kappa =  [0   0  ];
% loop1.phi1 =   [1.5 3.0];
% loop1.phi2 =   [2.0 3.5];

%exp5_5 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% outDir = 'exp5_5';
% param.R = 20;
% loop1.lambda = [0.1 0.1 0.1 0.1 ];
% loop1.kappa =  [0   0   0.1 0.1 ];
% loop1.phi1 =   [1.5 3.0 1.5 3.0 ]/10;
% loop1.phi2 =   [2.0 3.5 2.0 3.5 ]/10;

%exp5_6 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% outDir = 'exp5_6';
% param.R = 20;
% loop1.lambda = [0.1 0.1];
% loop1.kappa =  [0   0  ];
% loop1.phi1 =   [1.5 3.0]*2;
% loop1.phi2 =   [2.0 3.5]*2;

%exp5_7 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% opt.MaxMainIter = 20;
% outDir = 'exp5_7';
% param.R = 20;
% loop1.lambda = [0.1 0.1 0.1 0.1  0.1  0.1  0.1];
% loop1.kappa =  [0   0   0   0    0    0    0];
% loop1.phi1 =   [1.5 3.0 0.1 0.05 0.02 0.01 0]*2;
% loop1.phi2 =   [1.5 3.0 0.1 0.05 0.02 0.01 0]*2;

%exp6 rand uniform dict init,1max inner iters,param.Uinit = 'rand';
% Y = genSpotSeq2noisenorm;
% opt.MaxMainIter = 0;
% outDir = 'exp6';
% param.R = 20;
% loop1.lambda = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
% loop1.kappa =  [0   0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
% loop1.phi1 =   [0   0   0   0   0   0   0   0   0   0   0];
% loop1.phi2 =   [0   0   0   0   0   0   0   0   0   0   0];