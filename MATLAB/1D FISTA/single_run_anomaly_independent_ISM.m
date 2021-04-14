% Parameter selection
disp('Setup params')

% Parent directory
top_dir = 'D:\MMPAD_data_nr1\';
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\CHESS_data\';
% top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring4_zero';

% Output dirs
output_name = '_indep_ISM_mask1';
output_subdir = [dset_name,output_name];


% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)

num_ims = numel(dir(fullfile(dataset,'*.mat')));

% File Parameters
P.prefix = 'polar_vector';
P.baseFileName = 'indep_fit_%i_%i.mat';
P.dataset = dataset;

% Data/Dictionary Parameters
% Zero padding and mask

load(fullfile(dataset,[P.prefix,'_1.mat']));
N = numel(polar_vector(1:179))+180;
K = 20;
M = 30;
T = num_ims;
zPad = 0;
zMask = [];
% zMask = [1:zPad,(N+zPad+1):(N+2*zPad)];

P.dataScale = 1;
P.lambda_values = logspace(-4,1,M);
P.num_theta = N;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = linspace(0.5,50,P.num_var_t).^2;

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
P.params.plotProgress = 1;
P.params.verbose = 1;

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Load data
mPad = 90;

for j = 150
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = P.dataScale*sum(b_data.polar_image,1)';
    catch
        b = P.dataScale*polar_vector(1:179);
    end
    b = mirrorPad(b,mPad);
    B(:,j) = b;
end

%% Run coupled grid search
disp('Begin grid search')

% Init solution
x_init = zeros(size(A0ft_stack));

for i = 13
    P.set = i;
    P.params.lambda1 = P.lambda_values(i);
    for t = 19
        % Solve
        test_b = forceMaskToZero(B(:,t),zMask);
        [x_hat,obj,err,l1_norm] = convADMM_LASSO_Sherman_1D(A0ft_stack,test_b,x_init,P.params);  

        % Output data
        save_output(output_dir,x_hat,err,obj,l1_norm,P,t);
    end
end

%% Plot fit 
figure(1)
fit = Ax_ft_1D(A0ft_stack,x_hat);
hold on
plot(fit)
plot(test_b)
legend('fit','b')

%%
function save_output(output_dir,x_hat,err,obj,l1_norm,P,t)
    save(fullfile(output_dir,sprintf(P.baseFileName,P.set,t)),'x_hat',...
        'err','obj','l1_norm','P');
end


