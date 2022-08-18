% Parameter selection
disp('Setup params')
P.set = 1;

% Parent directory
top_dir = 'D:\MMPAD_data';
%     top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring1_zero_subset';

% Output dirs
output_name = '_indep_CG_TVphi1';
output_subdir = [dset_name,output_name];


% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)  

num_ims = numel(dir(fullfile(dataset,'*.mat')));

% File Parameters
P.prefix = 'mmpad_img';
P.baseFileName = 'indep_fit_%i.mat';
P.dataset = dataset;

% Data/Dictionary Parameters
% Zero padding and mask
zPad = [0,0];
zMask = [];
load(fullfile(dataset,[P.prefix,'_1.mat']));
polar_vector = sum(polar_image,1)';

N = size(polar_vector,1);
K = 20;
T = num_ims;

P.num_theta = N;
P.dtheta = 1;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = linspace(P.dtheta/2,50,P.num_var_t).^2;

% algorithm parameters

P.params.rho1 = 1;
P.params.lambda1 = 0.0001;
P.params.tau = 1.05;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 800;
P.params.conjGradIter = 50;
P.params.tolerance = 1e-6;
P.params.cgEpsilon = 1e-1;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end

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

for i = 1
    P.set = 2;
    
    % Init solution
    X_init = zeros(N,K,T);
    
    % Solve
    [X_hat,err,obj,l1_norm,tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack,B,X_init,P.params);  

    % Output data
    save_output(output_dir,X_hat,err,obj,l1_norm,tv_penalty,P);
end


function save_output(output_dir,X_hat,err,obj,l1_norm,tv_penalty,P)
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set)),'X_hat',...
        'err','obj','l1_norm','tv_penalty','P');
end


