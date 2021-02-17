% Parameter selection
disp('Setup params')

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'E:\MMPAD_data';
top_dir = '/cluster/shared/dbanco02';
noise_added = 0:0.03:0.30;
noise_thresh = [0.01,0.03:0.03:0.24,0.3,0.5];

for iii = 1:11
% Input dirs
dset_name = ['simulated_two_spot_1D_anomaly_small_',num2str(iii)];

% Indep dirs
indep_name = '_indep_ISM1';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
output_name = '_coupled_CG_TVphi2';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)  

% File Parameters
baseFileName = 'indep_fit_%i_%i.mat';
load(fullfile(indep_dir,sprintf(baseFileName,1,1)));
P.baseFileName = 'coupled_fit_%i.mat';

N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;

% Function
funcName = 'wrap_convADMM_LASSO_CG_TVphi_1D';

%% Fixed Parameters
% Zero padding and mask
zPad = [0,0];
zMask = [];

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Algorithm parameters
P.params.rho2 = 0.005;
P.params.lambda2 = 1;
P.params.tau = 1.1;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 800;
P.params.conjGradIter = 50;
P.params.tolerance = 1e-8;
P.params.cgEpsilon = 1e-1;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Lambda1 values: Use L-curve parameter selection
indep_data = load(fullfile(indep_dir,sprintf(baseFileName,1,1)));
lambda1_vals = indep_data.P.lambda_values;
M_lam1 = numel(lambda1_vals);
err_select = zeros(M_lam1,T);
l1_select = zeros(M_lam1,T);

for m = 1:M_lam1
    for t = 1:T
        load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
        b = P.dataScale*polar_vector(1:179);
        x_data = load(fullfile(indep_dir,sprintf(baseFileName,m,t)),'x_hat');
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_data.x_hat),129:133);
        err_select(m,t) = sum((fit(:)-b(:)).^2);
        l1_select(m,t) = sum(x_data.x_hat(:));
    end
end

%% L curve parameter selection for l1-norm term
select_indices = zeros(T,1);
for t = 1:T
    err_t = err_select(:,t);
    l1_t = l1_select(:,t);
    err_t = err_t;%/max(err_t);
    l1_t = l1_t;%/max(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist + (err_t == 0) )  );
end

for t = 1:T
    load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    b = P.dataScale*polar_vector(1:179);
    rel_err_t = err_select(:,t)/sum(b(:).^2);
    while rel_err_t(select_indices(t)) > noise_thresh(data_num)
        if select_indices(t) > 1
            select_indices(t) = select_indices(t) - 1;
        else
            select_indices(t) = find(rel_err_t == min(rel_err_t));
            break
        end
    end
end

P.params.lambda1 = lambda1_vals(select_indices);
P.params.lambda1_indices = select_indices;

% Lambda2 values
M = 30;
lambda2_vals = logspace(-4,1,M);
M = numel(lambda2_vals);
P.lambda2_values = lambda2_vals;

%% Parameters to vary

% Job directory
jobDir = fullfile('/cluster','home','dbanco02',['job_',output_subdir]);
mkdir(jobDir)

for k = 1:M
    P.params.lambda2 = lambda2_vals(k);
    P.set = k;
    varin = {dataset,P,output_dir};
    save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
end

slurm_write_bash(k-1,jobDir,'full_batch_script.sh',['1-',num2str(M)])
% slurm_write_matlab(k-1,jobDir,'parallel_FISTA','matlab_batch_script.sh')
end
