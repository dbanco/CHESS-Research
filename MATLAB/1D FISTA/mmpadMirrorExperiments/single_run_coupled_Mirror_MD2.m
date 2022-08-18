% Parameter selection
disp('Setup params')

% Parent directory
top_dir = 'D:\MMPAD_data_nr1';
% top_dir = 'E:\MMPAD_data';
% top_dir = '/cluster/shared/dbanco02';

ring_num = 1
% Input dirs
dset_name = ['ring',num2str(ring_num),'_zero'];

% Indep dirs
indep_name = '_indep_ISM_Mirror5';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
output_name = '_coupled_CG_TVphi_Mirror22';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
mkdir(output_dir)  

% File Parameters
indepName = 'indep_fit_%i_%i.mat';
load(fullfile(indep_dir,sprintf(indepName,30,1)));
P.baseFileName = 'coupled_fit_%i.mat';

N = P.num_theta;
K = P.num_var_t;
T = 10;

% Function
funcName = 'wrap_convADMM_LASSO_CG_TVphi_Mirror_1D';

% Fixed Parameters
% Zero padding and mask
zPad = [0,0];
zMask = [];

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

% Algorithm parameters
M = 30;
lambda2_vals = logspace(-5,1,M);
P.params.rho2 = 1e-5;
P.params.lambda2 = 5e-1;
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
P.params.plotProgress = 1;
P.params.verbose = 1;

% Load data
B = zeros(N,T);
for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = sum(b_data.polar_image,1);
        b(129:133) = (b(128) + b(134))/2;
    catch
        b = b_data.polar_vector;
    end

    % Mirror data
    nn = numel(b);
    pad1 = floor(nn/2);
    pad2 = ceil(nn/2);
    N = nn + pad1 + pad2;
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+nn)) = b;
    b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flip(b(1:pad1));    
    B(:,j) = b_mirror;

end

P.dataScale = 1;
P.dataScale = 1/mean(B(:));
B = B*P.dataScale;

%% Sparsity parameter selection
indep_data = load(fullfile(indep_dir,sprintf(indepName,1,1)));
lambda1_vals = indep_data.P.lambda_values;
M_lam1 = numel(lambda1_vals);
err_select = zeros(M_lam1,T);
l1_select = zeros(M_lam1,T);

rho1_select = zeros(M_lam1,T);
for m = 1:M_lam1
    fprintf('%i of %i \n',m,M_lam1)
    for t = 1:T
        x_data = load(fullfile(indep_dir,sprintf(indepName,m,t)),'x_hat','rho');
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_data.x_hat),129:133);
        err_select(m,t) = sum( (fit(:)-B(:,t)).^2 )/norm(B(:,t));
        l1_select(m,t) = sum(x_data.x_hat(:));
        rho1_select(m,t) = x_data.rho;
    end
end

% L curve parameter selection for l1-norm term
select_indices = zeros(T,1);
for t = 1:T
    err_t = err_select(:,t);
    l1_t = l1_select(:,t);
    err_t = err_t/max(err_t);
    l1_t = l1_t/max(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist + (err_t == 0) )  );
end

% Select minimum rho value
rho1 = max(rho1_select(:));
P.params.lambda1 = zeros(546,1);
for t = 1:T
    rho1 = min(rho1_select(select_indices(t),t),rho1);
    P.params.lambda1(t) = lambda1_vals(select_indices(t))*1e-8;
end
P.params.rho1 = rho1;

P.params.lambda1(:) = mean(P.params.lambda1(t));
%% Run coupled algorithm
% Init solution
P.params.rho = 1e-2;
X_init = zeros([size(A0ft_stack),T]);
for t = 1:T
    i_data = load(fullfile(indep_dir,sprintf(indepName,select_indices(t),t)));
    X_init(:,:,t) = i_data.x_hat;    
end

%Solve
[X_hat,err,obj,l1_norm,tv_penalty] = convADMM_LASSO_CG_TVphi2_1D(A0ft_stack,B,X_init,P.params); 


%% Plot AWMV
P.var_theta = [linspace(0.5,100,30)].^2;
awmv = zeros(T,1);
var_signal = squeeze(sum(X_hat,1));
var_sum = squeeze(sum(var_signal,1));
for t = 1:T
    awmv(t) = sum(sqrt(P.var_theta(:)).*var_signal(:,t))/var_sum(t);
end

figure(1)
plot(awmv)

%% Plot fit 
Bnorms = zeros(T,1);
B_normed = zeros(size(B));
for j = 1:T
    Bnorms(j) = norm(B(:,j));
    B_normed(:,j) = B(:,j)/Bnorms(j);
end
fit = Ax_ft_1D_Time(A0ft_stack,X_hat,Bnorms);

figure(2)
subplot(2,1,1)
imagesc(B_normed)
colorbar()
subplot(2,1,2)
imagesc(fit)
colorbar()

figure(3)
im_num = 1;
hold on
plot(B_normed(:,im_num))
plot(fit(:,im_num))

%% Plot fits for L-curve selected parameter
fits_fig = figure(9);
[ha2, ~] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
for t = 1:10   
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(B_normed(:,t))
    plot(fit(:,t))
    rel_err = sum(( fit(:,t) - B_normed(:,t) ).^2);
    legend(sprintf('%i \n %0.2f',sum(x_hat(:)>0),rel_err),'location','northeast')
    im_ind = im_ind + 1;
end

