% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_polar_fit_alpha_beta');
mkdir(outputdir)
% Function
funcName = 'wrap_FISTA_Circulant';

jobDir = fullfile(datadir,'job_al7075_311_alpha_beta');
mkdir(jobDir)

%% Fixed Parameters
P.img = 35;
P.load_step = 0;

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  3,    P.num_var_r).^2;

% basis weighting
P.weight = 1;

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
alphaps = [1 5 10 20 50 100 200 500 1000 5000 10000 100000 1000000];
betaps = [0.001 0.01 0.5 0.1 1 5 10 100 1000]*P.dtheta*P.drad;

k = 0;
for alphap = alphaps
    P.alphap = alphap;
    for betap = betaps
        P.betap = betap;
        P.task = k;
        varin = {dataset,P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

slurm_write_bash(k-1,jobDir)
