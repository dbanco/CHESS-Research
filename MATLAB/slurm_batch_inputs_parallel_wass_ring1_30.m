% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'mmpad_polar');
ringName = 'ring1_zero';
ring_num  = 1;
prefix = 'mmpad_img';

% Output directories
dirA = fullfile('/cluster','shared','dbanco02','wass_parallel_a_30');
dirB = fullfile('/cluster','shared','dbanco02','wass_parallel_b_30');
mkdir(dirA)
mkdir(dirB)

% Functions
funcName1 = 'wrap_FISTA_Circulant';
funcName2 = 'wrap_wass_reg_FISTA_Circulant';
jobDir1 = fullfile(datadir,['job_wass_parallel_init_30']);
jobDir2 = fullfile(datadir,['job_wass_parallel_ab_30']);
jobDir3 = fullfile(datadir,['job_wass_parallel_ba_30']);
mkdir(jobDir1)
mkdir(jobDir2)
mkdir(jobDir3)

%% Universal Parameters
% Ring sampling parameters
load(fullfile(dataset,ringName,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [31,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 15;
P.num_var_r = 8;
P.var_theta   = linspace(P.dtheta/2,  32,P.num_var_t).^2;
P.var_rad = linspace(P.drad/2,6,P.num_var_r).^2;

% Zero padding and mask
maskCols = 129:133;
zPad = [0,0];
zMask = zeros(size(zeroPad(polar_image,zPad)));
zMask(:,maskCols) = 1;
zMask = onePad(zMask,zPad);
[r,c] = find(zMask==1);
zMask = [r,c];


%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 10;
params.t_k = 1;
params.lambda = 0.01;
params.wLam = 30;
params.gamma = 0.008;
params.beta = 1.2;
params.maxIter = 200;
params.maxIterReg = 5;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Initialization parameters
P1 = P;
P1.params.tolerance = 1e-8;

%% Regularization parameters
P2 = P;

%% Parameters to vary
img_nums = 20:50;

imageDir = fullfile(dataset,ringName,prefix);
% Initialization jobs
k = 1;
for img = img_nums
    P1.img = img;
    P1.index = k;
    P1.set = ring_num;
    varin = {imageDir,P1,dirA};
    funcName = funcName1;
    save(fullfile(jobDir1,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

% Regularization jobs for dirA->dirB
k = 1;
for img = img_nums
    P2.img = img;
    P2.index = k;
    P2.set = ring_num;
    varin = {dirA,P2,dirB};
    funcName = funcName2;
    save(fullfile(jobDir2,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

% Regularization jobs for dirB->dirA
k = 1;
for img = img_nums
    P2.img = img;
    P2.index = k;
    P2.set = ring_num;
    varin = {dirB,P2,dirA};
    funcName = funcName2;
    save(fullfile(jobDir3,['varin_',num2str(k),'.mat']),'varin','funcName')
    k = k + 1;
end

% Init script
slurm_write_matlab(numel(img_nums),jobDir1,'parallel_reg_wass_FISTA_30','batch_script.sh')

