%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\Simulations';
top_dir = '/cluster/home/dbanco02/data/Simulations';
mkdir(top_dir)

% Simulation name
sim = 'linear';
sim_name = [sim,'PoissonNoise2'];
file_name = 'lineSearchNoise';

% Define poisson dataset
N = 101;
T = 30;
[P,K,M,MM] = definePoissonP(N,T);
levels = 0.05:0.05:0.3;
[alpha_vals,theta_stds1] = genPoissonScales(N,T,levels,sim);
P.theta_stds = theta_stds1;

%% Run algorithm
% Independent
alg_name = 'IndepISM';
indep_dir  = fullfile(top_dir,sim_name,alg_name);
f_name =  [file_name,'_',num2str(nn),'.mat'];
       
% Construct dictionary
A0 = unshifted_basis_vector_ft_stack_zpad(P);

% Generate data
if strcmp(sim,'linear')
    [B,~,theta_stds1,~] = genLinearPoisson(N,T,alpha_vals(nn));
elseif strcmp(sim,'anomaly')
    [B,~,theta_stds1,~] = genAnomalyPoisson(N,T,alpha_vals(nn));
end

% Independent Solution
x_init = zeros(N,K);
X_indep = zeros(N,K,M,T);
i = 48;
P.set = i;
P.params.lambda1 = 0.1*P.lambda_values(i);
t = 22;
% Solve
P.params.plotProgress = 0;
P.params.rho1 = 1;
[x_hat] = convADMM_LASSO_Sherman_1D(A0,B(:,t),x_init,P.params);  
fit = Ax_ft_1D(A0,x_hat);

opt.Verbose = 1;
opt.MaxMainIter = 50;
opt.Y0 = x_init;
opt.rho = P.params.rho1;
opt.AutoRho = 1;
opt.NonNegCoef = 1;
[x_hat2, optinf] = cbpdn(ifft(A0), B(:,t), P.params.lambda1, opt);

fit2 = Ax_ft_1D(A0,x_hat2);
figure(1)
hold on
plot(B(:,t))
plot(fit)
plot(fit2)
legend('data','fit','fit2')
