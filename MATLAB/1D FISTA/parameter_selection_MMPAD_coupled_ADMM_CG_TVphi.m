% Parameter selection
disp('Setup params')

% Parent directory
top_dir = 'D:\MMPAD_data';
%     top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring1_zero_subset';

% Indep dirs
indep_name = '_indep_ISM3';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
output_name = '_coupled_CG_TVphi1';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)  

% File Parameters
P.prefix = 'mmpad_img';
P.baseFileName = 'coupled_fit_%i.mat';
P.dataset = dataset;

% Data/Dictionary Parameters
% Zero padding and mask
zPad = [0,0];
zMask = [];
load(fullfile(dataset,[P.prefix,'_1.mat']));
polar_vector = sum(polar_image,1)';

N = size(polar_vector,1);
K = 20;
M = 30;
T = 10;

P.num_theta = N;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = linspace(P.dtheta/2,50,P.num_var_t).^2;

% algorithm parameters
P.params.rho1 = 1;
P.params.lambda1 = [];
P.params.rho2 = 1;
P.params.lambda2 = 1;
P.params.tau = 1.05;
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
baseFileName = 'indep_fit_%i_%i.mat';
indep_data = load(fullfile(indep_dir,sprintf(baseFileName,1,1)));
lambda1_vals = indep_data.P.lambda_values;
select_indices = zeros(T,1);
M_lam1 = numel(lambda1_vals);
err_select = zeros(M_lam1,T);
l1_select = zeros(M_lam1,T);

for m = 1:M_lam1
    for t = 1:T
        load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
        b = sum(polar_image,1);
        e_data = load(fullfile(indep_dir,sprintf(baseFileName,m,t)),'err','x_hat');
        err_select(m,t) = e_data.err(end)/norm(b).^2;
        l1_select(m,t) = sum(e_data.x_hat(:))/norm(b);
        
    end
end

for t = 1:T
    err_t = err_select(:,t);
    l1_t = l1_select(:,t);
    sq_origin_dist = l1_t.^2 + err_t.^2;
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist) );
end

P.params.lambda1 = select_indices;

% Lambda2 values
lambda2_vals = logspace(1,5,M);
M = numel(lambda2_vals);
P.lambda2_values = lambda2_vals;

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);

% Load data
B = zeros(N,T);
for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = sum(b_data.polar_image,1);
    catch
        b = b_data.polar_vector;
    end
    B(:,j) = b';
end

%% Run coupled grid search
disp('Begin grid search')
for i = 1:M
    P.params.lambda2 = lambda2_vals(i);
    P.set = i;
    
    % Init solution
    X_init = zeros(N,K,T);
    
    % Solve
    [X_hat,err,obj,l1_norm,tv_penalty] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,B,X_init,P.params);  

    % Output data
    save_output(output_dir,X_hat,err,obj,l1_norm,tv_penalty,P);
end


function save_output(output_dir,X_hat,err,obj,l1_norm,tv_penalty,P)
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set)),'X_hat',...
        'err','obj','l1_norm','tv_penalty','P');
end


