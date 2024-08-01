% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'al7075_311_polar');

% Output directory
outputdir = fullfile('/cluster','shared','dbanco02','al7075_311_polar_fit');
mkdir(outputdir)
% Function
funcName = 'wrap_Fls
ISTA_Circulant';

jobDir = fullfile(datadir,'job_al7075_311');
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
alphaps = logspace(log10(0.005),log10(1),20);
betaps = logspace(log10(0.005),log10(10),20);

k = 0;
for alphap = alphaps
    for betap = betaps
        P.alphap = alphap;
        P.betap = betap;
        P.task = k;
        varin = {dataset,P,outputdir};
        save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

slurm_write_bash(k-1,jobDir)
