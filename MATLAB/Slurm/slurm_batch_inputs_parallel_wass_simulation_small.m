% Data directory
datadir = fullfile('/cluster','home','dbanco02');

% Ring dataset
dataset = fullfile(datadir,'simulated_data_small');
prefix = 'polar_image';

% Functions
funcName2 = 'wrap_wass_reg_FISTA_Circulant';


%% Universal Parameters
% Ring sampling parameters
load(fullfile(dataset,[prefix,'_1.mat']));
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;
P.sampleDims = [20,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  5,P.num_var_r).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];


%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-10;
params.L = 1000;
params.t_k = 1;
params.lambda = 0.0359;
params.wLam = 25;
params.gamma = 0.01;
params.beta = 1.2;
params.maxIter = 200;
params.maxIterReg = 5;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

%% Parameters to vary
img_nums = 1:20;
ring_num = 1;
gamma_vals = [0.00001 0.00002 0.00005 0.0001 0.0002 0.0005  0.001  0.002,...
                0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5];
imageDir = fullfile(dataset,prefix);

% Regularization jobs for dirA->dirB
k = 1;
for i = 1:numel(gamma_vals)
    
    inputdir = fullfile('/cluster','shared','dbanco02',...
                        ['wass_small_fit_a_',num2str(i)]);
    % Output directory
    outputdir = fullfile('/cluster','shared','dbanco02',...
                        ['wass_small_fit_b_',num2str(i)]);
    mkdir(outputdir)
    % Job directory
    jobDir2 = fullfile(datadir,['job_wass_parallel_small_',num2str(i),'_ab']);
    mkdir(jobDir2)
    slurm_write_matlab(numel(img_nums),jobDir2,'parallel_small_wass_FISTA','batch_script.sh',sprintf('wass_parallel_small_%i',i))
    for img = img_nums
        P.img = img;
        P.index = img;
        P.set = i;
        P.params.gamma = gamma_vals(i);
        varin = {inputdir,P,outputdir};
        funcName = funcName2;
        save(fullfile(jobDir2,['varin_',num2str(img),'.mat']),'varin','funcName')
        k = k + 1;
    end
end

% Regularization jobs for dirB->dirA
k = 1;
for i = 1:numel(gamma_vals)
    inputdir = fullfile('/cluster','shared','dbanco02',...
                        ['wass_small_fit_b_',num2str(i)]);
    % Output directory
    outputdir = fullfile('/cluster','shared','dbanco02',...
                        ['wass_small_fit_a_',num2str(i)]);
    mkdir(outputdir)
    % Job directory
	jobDir3 = fullfile(datadir,['job_wass_parallel_small_',num2str(i),'_ba']);
    mkdir(jobDir3)
    for img = img_nums
        P.img = img;
        P.index = img;
        P.set = i;
        P.params.gamma = gamma_vals(i);
        varin = {inputdir,P,outputdir};
        funcName = funcName2;
        save(fullfile(jobDir3,['varin_',num2str(img),'.mat']),'varin','funcName')
        k = k + 1;
    end
end
% Init script


