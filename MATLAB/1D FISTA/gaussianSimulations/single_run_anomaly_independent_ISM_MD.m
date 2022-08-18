close all

% Parameter selection
disp('Setup params')

% Parent directory
top_dir = 'D:\MMPAD_data_nr1\';
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\CHESS_data\';
% top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring1_zero';

% Output dirs
output_name = '_indep_ISM_mask11';
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

load(fullfile(dataset,[P.prefix,'_1.mat']));
polar_vector = sum(polar_image,1);
n = numel(polar_vector);
pad1 = floor(n/2);
pad2 = ceil(n/2);
N = n + pad1 + pad2;
K = 30;
M = 30;
T = num_ims;

zPad = 0;
zMask = [];
% zMask = [1:zPad,(n+zPad+1):(n+2*zPad)];

P.lambda_values = logspace(-4,1,M);
P.num_theta = N;
P.sampleDims = [T,1];
P.num_ims = T;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = linspace(0.5,100,P.num_var_t).^2;

% algorithm parameters
P.params.rho1 = 0.00001;
% P.params.lambda1 = 0.0001;
P.params.tau = 1.05;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 800;
P.params.tolerance = 1e-10;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Load data
B = zeros(n ,T);
B_m = zeros(N ,T);
Bnorms = zeros(T ,1);

for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    % Reduce image to vector if needed
    try
        b = sum(b_data.polar_image,1)'*1e-4;
        b(129:133) = (b(128) + b(134))/2;
    catch
        b = polar_vector;
    end
    nn = numel(b);
    B(:,j) = b;
    
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+n)) = b;
    b_mirror((1+N-pad2):N) = flipud(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flipud(b(1:pad1));
    B_m(:,j) = b_mirror;
end

% P.dataScale = 1;
% P.dataScale = 1/mean(B(:));
% P.dataScale = 1e-4;
% B = B*P.dataScale;

%% Run coupled grid search
% disp('Begin grid search')
% 
% % Init solution
% x_init = zeros(size(A0ft_stack));
% load('lam150')
% 
% t = 150;
% test_b = forceMaskToZero(B(:,t),zMask);
% 
% % Solve
% P.params.lambda1 = 2e-4;
% b1 = test_b;
% [x_hat1,~,~,~] = convADMM_LASSO_MaskDecouple_1D(A0ft_stack/norm(b1),b1/norm(b1),x_init,P.params);
% P.params.lambda1 = 8e-2;
% b2 = test_b;
% [x_hat2,~,~,~] = convADMM_LASSO_MaskDecouple_1D(A0ft_stack,b2,x_init,P.params);

%% Run coupled grid search
disp('Begin grid search')

%Init solution

load('lam150')

t = 101;
test_b = B(:,t);

%Solve
lam1 = 1;
P.params.lambda1 = lam1;

b1 = B(:,t);
P.num_theta = n;
A0ft_stack1 = unshifted_basis_vector_ft_stack_zpad(P);
x_init = zeros(size(A0ft_stack1));
[x_hat1,~,~,~] = convADMM_LASSO_Sherman_1D(A0ft_stack1/norm(b1),b1/norm(b1),x_init,P.params);

b2 = B_m(:,t);
P.num_theta = N;
P.params.maxIter = 800;
lam2 = lam1;
P.params.lambda1 = lam2;
A0ft_stack2 = unshifted_basis_vector_ft_stack_zpad(P);
x_init = zeros(size(A0ft_stack2));
[x_hat2,~,~,~] = convADMM_LASSO_Sherman_1D(A0ft_stack2/norm(b1),b2/norm(b1),x_init,P.params);

% P.params.lambda1 = 2e-4*norm(b1).^2;
% b2 = test_b;
% [x_hat2,~,~,~] = convADMM_LASSO_Sherman_1D(A0ft_stack,b2,x_init,P.params);

%% Plot fit 

% x_hat1 = zeros(size(x_hat1));
% x_hat1(end) = 50;
fit1 = Ax_ft_1D(A0ft_stack1,x_hat1);
fit2 = Ax_ft_1D(A0ft_stack2,x_hat2);
fprintf('              Unaltered,   Mirrored,     Mirrored-cropped\n')

fprintf('lambda:       %2.2e,    %2.2e,     %2.2e \n',....
                    lam1,...
                    lam2,...
                    lam2 )

fprintf('Error:        %2.2e,,   %2.2e,     %2.2e\n',...
                    norm(fit1-b1)^2/norm(b1)^2,...
                    norm(fit2-b2)^2/norm(b2)^2,...
                    norm(fit2(pad1+1:pad1+n)-b2(pad1+1:pad1+n))^2/norm(b1)^2) 
                
fprintf('||x||_1:      %2.2e,    %2.2e,     %2.2e\n',....
                    sum(abs(x_hat1),'all'),...
                    sum(abs(x_hat2),'all'),...
                    sum(abs(x_hat2( (pad1+1:pad1+n),:)),'all')  )
fprintf('lam*||x||_1:  %2.2e,    %2.2e,     %2.2e\n',...
                    lam1*sum(abs(x_hat1),'all'),...
                    lam2*sum(abs(x_hat2),'all'),...
                    lam2*sum(abs(x_hat2( (pad1+1:pad1+n),:)),'all')  )
figure(1)
subplot(3,1,1)
hold on
plot(b1)
plot(fit1)
legend('data','fit1')

subplot(3,1,2)
hold on
plot(b2(pad1+1:pad1+n))
plot(fit2(pad1+1:pad1+n))
legend('data','fit2 (cropped)')

subplot(3,1,3)
hold on
plot(b2 )
plot(fit2 )
legend('data','fit2')



%%


