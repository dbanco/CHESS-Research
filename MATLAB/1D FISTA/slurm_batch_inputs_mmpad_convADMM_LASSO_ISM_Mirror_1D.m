% Parameter selection
disp('Setup params')

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'E:\MMPAD_data';
top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring4_zero';

% Indep dirs
indep_name = '_indep_ISM_Mask1';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);
mkdir(indep_dir)

% Setup directories
dataset =  fullfile(top_dir,dset_name);

num_ims = numel(dir(fullfile(dataset,'*.mat')));

% File Parameters
P.prefix = 'mmpad_img';
P.baseFileName = 'indep_fit_%i_%i.mat';
P.dataset = dataset;

% Data/Dictionary Parameters
% Zero padding and mask
load(fullfile(dataset,[P.prefix,'_1.mat']));
polar_vector = sum(polar_image,1);
n = numel(polar_vector);
zPad = 0;
% zMask = [1:zPad,(zPad+n+1):(2*zPad+n)];
zMask = [];

K = 20;
M = 50;
T = num_ims;
N = n + floor(n/2) + ceil(n/2);

P.dataScale = 1;
P.lambda_values = logspace(-5,-1,M);
P.num_theta = N;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = [linspace(0.5,50,P.num_var_t)].^2;

% algorithm parameters
P.params.rho1 = 0.1;
% P.params.lambda1 = 0.0001;
P.params.tau = 1.05;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 1000;
P.params.tolerance = 1e-12;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Function
funcName = 'wrap_convADMM_LASSO_ISM_Mirror_1D';

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_',indep_subdir]);
mkdir(jobDir)

for k = 1:M
    P.set = k;
    P.params.lambda1 = P.lambda_values(k);
    
    varin = {dataset,P,indep_dir};
    save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
end

slurm_write_bash(M,jobDir,'full_batch_script.sh',['1-',num2str(M)])
% slurm_write_matlab(k-1,jobDir,'parallel_FISTA','matlab_batch_script.sh')
