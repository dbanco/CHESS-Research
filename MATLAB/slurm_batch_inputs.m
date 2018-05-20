% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar_reduced');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_norm2_2');

mkdir(outputdir)
% Function
funcName = 'wrap_FISTA_Circulant';

jobDir = fullfile('/cluster','home','dbanco02','job_al7075_311_norm2_2');
mkdir(jobDir)

%% Fixed Parameters

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;
P.sampleDims = [37,5];

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,    P.num_var_r).^2;

% basis weighting
P.weight = 1;
P.alphap = 1;
P.betap = 1;

% fista params
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 500;
params.isNonnegative = 1;
P.params = params;

%% Parameters to vary
img_nums = 0:184;
load_steps = 0:4;

k = 0;
for load_step = load_steps
    for img = img_nums
        P.img = img;
        P.load_step = load_step;

        varin = {dataset,P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

slurm_write_bash(k-1,jobDir,'full_batch_script.sh','0-924')