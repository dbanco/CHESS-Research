% Parameter selection
disp('Setup params')

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'E:\MMPAD_data';
top_dir = '/cluster/shared/dbanco02/data/MMPAD_omega';
om_dir = {'omega2','omega3','omega4','omega5'};
r_dir = {'ring1','ring2','ring3','ring4'};


for o = 4
for ring_num = 1
% Input dirs
dset_name = ['ring',num2str(ring_num)];

% Indep dirs
indep_name = '_indep_ISM_Mirror5';
indep_subdir = [dset_name,om_dir{o},indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
output_name = '_coupled_CG_TVphi_Mirror5';
output_subdir = [dset_name,om_dir{o},output_name];

% Setup directories
dataset =  fullfile(top_dir,om_dir{o},dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)  

% File Parameters
baseFileName = 'indep_fit_%i_%i.mat';
load(fullfile(indep_dir,sprintf(baseFileName,1,1)));
P.baseFileName = 'coupled_fit_%i.mat';

P.indepDir = indep_dir;
P.indepName = baseFileName;
P.indepInit = 1;

N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;
T = 546;
P.num_ims = 546;

% Function
funcName = 'wrap_convADMM_LASSO_CG_TVphi_Mirror_1D';

%% Fixed Parameters
% Zero padding and mask
zPad = [0,0];
zMask = [];

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Algorithm parameters
P.params.rho2 = 1e-4;
P.params.lambda2 = 1;
P.params.tau = 1.1;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 50;
P.params.conjGradIter = 50;
P.params.tolerance = 1e-8;
P.params.cgEpsilon = 1e-3;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Load data
B = zeros(N,T);
for t = 1:T
    load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    b = P.dataScale*sum(polar_image,2);
    % Mirror data
    nn = numel(b);
    pad1 = floor(nn/2);
    pad2 = ceil(nn/2);
    N = nn + pad1 + pad2;
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+nn)) = b;
    b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flip(b(1:pad1));
    B(:,t) = b_mirror;
end

% Lambda1 values: Use L-curve parameter selection
indep_data = load(fullfile(indep_dir,sprintf(baseFileName,1,1)));
lambda1_vals = indep_data.P.lambda_values;
M_lam1 = numel(lambda1_vals);
err_select = zeros(M_lam1,T);
l1_select = zeros(M_lam1,T);

rho1_select = zeros(M_lam1,T);
for m = 1:M_lam1
    for t = 1:T
        x_data = load(fullfile(indep_dir,sprintf(baseFileName,m,t)),'x_hat','rho');
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_data.x_hat),129:133);
        err_select(m,t) = sum( (fit(:)-B(:,t)).^2 );
        l1_select(m,t) = sum(x_data.x_hat(:));
        rho1_select(m,t) = x_data.rho;
    end
end

%% L curve parameter selection for l1-norm term
select_indices = zeros(T,1);
for t = 1:T
    err_t = err_select(:,t);
    l1_t = l1_select(:,t);
    err_t = log(err_t);
    l1_t = log(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist )  );
end

% for t = 1:T
%     load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
%     b = P.dataScale*sum(polar_image,1);
%     rel_err_t = err_select(:,t)/sum(b(:).^2);
%     while rel_err_t(select_indices(t)) > 0.02
%         if select_indices(t) > 1
%             select_indices(t) = select_indices(t) - 1;
%         else
%             select_indices(t) = find(rel_err_t == min(rel_err_t));
%             break
%         end
%     end
% end

% P.params.lambda1 = lambda1_vals(select_indices);
P.params.lambda1 = ones(T,1)*0.4546e-3;
P.params.lambda1_indices = select_indices;

% Select minimum rho value
% rho1 = max(rho1_select(:));
% for t = 1:T
%     rho1 = min(rho1_select(select_indices(t),t),rho1);
% end
% P.params.rho1 = rho1;

% Lambda2 values
M = 30;
lambda2_vals = logspace(-5,1,M);
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
end