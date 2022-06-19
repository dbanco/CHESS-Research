  % cbpdndl_sep_script
% Load MMPAD subset
data_dir = 'D:\MMPAD_data_nr1\ring1_zero';
t = 1:4:200;
T = numel(t);
for i = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t(i)),'.mat']))
    if t == 1
        [N,M] = size(polar_image(:,134:end));
        Y = zeros(N,M,T);
    end
    Y(:,:,i) = polar_image(:,134:end)/norm(polar_image(:,134:end),'fro');
end


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

% outDir = 'exp1';
% loop1.lambda = 0.01:0.01:0.2;

outDir = 'exp4_group';
% loop1.lambda = linspace(0.125,0.24,10);
loop1.lambda = 0.24*ones(11,1);
loop1.kappa = [0.001 0.005 0.01 0.05 1 0 10 100 1000 10000 100000];

saveOuts = 1;
topDir = ['C:\Users\dpqb1\Desktop\cbpdndl_sep_mmpad\',outDir];
mkdir(topDir)

for ii = 9:11
    ending = ['_out_',num2str(ii)];
    kappa = loop1.kappa(ii)
    lambda = loop1.lambda(ii);
    opt.rho = (50*lambda + 0.5);

    %% Solve
    [D, X, Dx, Dy, optinf] = cbpdndl_sep_group_smooth(D0, D0h, D0v, Y, lambda, kappa, opt);
%     load(fullfile(topDir,['variables',ending,'.mat']))
    saveOutputsCBPDNDL_SepSmooth(topDir,ending,Y,D, X, Dx, Dy, optinf,lambda)
    
end