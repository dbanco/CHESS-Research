% Parameter selection
disp('Setup params')

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'E:\MMPAD_omega';
% top_dir = '/cluster/home/dbanco02/data/MMPAD_omega';
om_dir = {'omega2','omega3','omega4','omega5'};
r_dir = {'ring1','ring2','ring3','ring4'};

awmv_all_indep = zeros(546,50,4);
all_select_indices = zeros(546,4);
all_select_indices_err = zeros(546,4); 

for o = 4
for ring_num = 2:4
    fprintf('%i, ',ring_num)
% Input dirs
dset_name = ['ring',num2str(ring_num)];

% Indep dirs
indep_name = '_indep_ISM_Mirror';
indep_subdir = [dset_name,om_dir{o},indep_name];
indep_dir = fullfile(top_dir,'indep',indep_subdir);

% Setup directories
dataset =  fullfile(top_dir,om_dir{o},dset_name);

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
P.params.rho2 = 0.0001;
P.params.lambda2 = 1;
P.params.tau = 1.1;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 1000;
P.params.conjGradIter = 50;
P.params.tolerance = 1e-8;
P.params.cgEpsilon = 1e-5;
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
rel_err_select = zeros(M_lam1,T);

for m = 1:M_lam1
    for t = 1:T
        x_data = load(fullfile(indep_dir,sprintf(baseFileName,m,t)),'x_hat','rho');
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_data.x_hat),129:133);
        err_select(m,t) = sum( (fit(:)-B(:,t)).^2 );
        rel_err_select(m,t) = norm(fit(:)-B(:,t))/norm(B(:,t));
        l1_select(m,t) = sum(x_data.x_hat(:));
        var_signal = sum(x_data.x_hat,1);
        var_total = sum(var_signal(:));
        awmv_all_indep(t,m,ring_num) = sum(var_signal.*sqrt(P.var_theta))/var_total;
        
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
all_select_indices(:,ring_num) = select_indices;


for t = 1:T
    rel_err_t = rel_err_select(:,select_indices(t));
    while rel_err_t > 0.02
        if select_indices(t) > 1
            select_indices(t) = select_indices(t) - 1;
        else
            select_indices(t) = find(rel_err_t == min(rel_err_t));
            break
        end
    end
end
all_select_indices_err(:,ring_num) = select_indices;
end
end

save('om5_awmv_indices_indep.mat','all_select_indices','all_select_indices_err','awmv_all_indep','err_select','rel_err_select','l1_select')
