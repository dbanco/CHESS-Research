% Parameter selection
disp('Setup params')

% Parent directory
top_dir = 'E:\MMPAD_data';
% top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring3_zero';

% Output dirs
output_name = '_indep_ISM4';
output_subdir = [dset_name,output_name];


% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)

num_ims = numel(dir(fullfile(dataset,'*.mat')));

% File Parameters
P.prefix = 'mmpad_img';
P.baseFileName = 'indep_fit_%i_%i.mat';
P.dataset = dataset;

% Data/Dictionary Parameters
% Zero padding and mask
zPad = 0;
zMask = [];
load(fullfile(dataset,[P.prefix,'_1.mat']));
polar_vector = sum(polar_image,1)';

N = size(polar_vector,1);
K = 20;
M = 30;
T = num_ims;

P.dataScale = 1e-5;
P.lambda_values = logspace(-5,1,M);
P.num_theta = N;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = logspace(-2,2,P.num_var_t).^2;

% algorithm parameters
P.params.rho1 = 1;
% P.params.lambda1 = 0.0001;
P.params.tau = 1.05;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 800;
P.params.tolerance = 1e-8;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Load data
B = zeros(N,T);
for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = P.dataScale*sum(b_data.polar_image,1);
    catch
        b = P.dataScale*b_data.polar_vector;
    end
    B(:,j) = b';
    B(129:133,j) = (b(128) + b(134))/2;
end

%% Run coupled grid search
disp('Begin grid search')

% Init solution
x_init = zeros(N,K);

for i = 1:M
    P.set = i;
    P.params.lambda1 = P.lambda_values(i);
    for t = 1:T
        % Solve
        [x_hat,obj,err,l1_norm] = convADMM_LASSO_Sherman_1D(A0ft_stack,B(:,t),x_init,P.params);  

        % Output data
        save_output(output_dir,x_hat,err,obj,l1_norm,P,t);
    end
end


function save_output(output_dir,x_hat,err,obj,l1_norm,P,t)
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set,t)),'x_hat',...
        'err','obj','l1_norm','P');
end


