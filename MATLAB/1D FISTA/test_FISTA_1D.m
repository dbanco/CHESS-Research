% Test 1D FISTA
inDir = 'D:\CHESS_data\simulated_two_spot_1D\polar_vector';
outDir = 'D:\CHESS_data\simulated_two_spot_1D_fit_single';
mkdir(outDir)

P.dtheta = 1;
P.sampleDims = [71,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 120;
P.var_theta   = linspace(P.dtheta/2,  50,P.num_var_t).^2;

% Zero padding and mask
zPad = [];
zMask = [];


% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 8000;
params.t_k = 1;
params.lambda = 0.06;
params.beta = 1.005;
params.maxIter = 3000;
params.maxIterReg = 50;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

for i = 1
    % Ring sampling parameters
    load([inDir,'_',num2str(i),'.mat']);
    P.img = i;
    P.set = 1;
    P.num_theta= size(polar_vector,2);

    wrap_FISTA_Circulant_1D(inDir,P,outDir)
end