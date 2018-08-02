% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar_synth2');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','synth2_norm2_1');
mkdir(outputdir)

% Function
funcName1 = 'wrap_norm2_FISTA_Circulant';
jobDir1 = fullfile(datadir,'job_synth2_norm2_1');
mkdir(jobDir1)

%% Fixed Parameters
% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;
P.sampleDims = [5,5];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  6,    P.num_var_r).^2;
% Zero padding
params.zeroPad = [10,0];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 1;
params.lambda = 0.1;
params.beta = 1.1;
params.maxIter = 500;
params.isNonnegative = 1;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary
k = 0;
for img = 0:24
    for load_step = 0
        P.img = img;
        P.load_step = load_step;
        varin = {dataset,P,outputdir};
        funcName = funcName1;
        save(fullfile(jobDir1,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

% Init script
slurm_write_bash(k-1,jobDir1,'fista_batch_script.sh','0-24')