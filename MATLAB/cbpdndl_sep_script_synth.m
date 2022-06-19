% cbpdndl_sep_script
%Generate toy data
Y = genSpotSeq2normPlusNoise;
[N,M,T] = size(Y);

% Set up cbpdndl parameters

opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.sigma = (size(Y,3));
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.LinSolve = 'CG';
opt.CGTol = 1e-3;
opt.CGTolAuto = 0;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;

% Construct initial dictionary
K = 8;
D0v = zeros(24,1,K);
D0v(:,:,:) = ones(24,1,K);

D0h = zeros(1,48,K);
D0h(:,:,:) = ones(1,48,K);

for c=1:size(D0h,3)
   D0(:,:,c) = D0v(:,:,c)*D0h(:,:,c);
end

% outDir = 'exp10';
% loop1.lambda = [0.01 0.02 0.03 0.04];

% Y now has additive noise w/ std = 1/100
outDir = 'exp11';
loop1.lambda = [0.01 0.02 0.03 0.04];

saveOuts = 1;
topDir = ['C:\Users\dpqb1\Desktop\cbpdndl_sep\',outDir];
mkdir(topDir)

for ii = 1:4
    ending = ['_out_',num2str(ii)];
    
    lambda = loop1.lambda(ii);
    opt.rho = (50*lambda + 0.5);

    %% Solve
    [D, X, Dx, Dy, optinf] = cbpdndl_sep2(D0, D0h, D0v, Y, lambda, opt);
%     load(fullfile(topDir,['variables',ending,'.mat']))
    saveOutputsCBPDNDL_Sep(topDir,ending,Y,D, X, Dx, Dy, optinf,lambda)
    
end