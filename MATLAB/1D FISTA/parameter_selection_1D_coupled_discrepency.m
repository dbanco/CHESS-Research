%% Parameter selection
clear all
close all
P.set = 1;
datadir = 'D:\CHESS_data\';
dataset = ['D:\CHESS_data\simulated_two_spot_1D_noise2_',dataset_num,'\'];
param_dir = 'D:\CHESS_data\param_search_1D\';
num_ims = 10;

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load([dataset,prefix,'_1.mat']);
P.num_theta = size(polar_vector,1);
P.dtheta = 1;
P.sampleDims = [num_ims,1];

% Basis function variance parameters
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = 15;
P.var_theta = linspace(P.dtheta/2,30,P.num_var_t).^2;

% Zero padding and mask\
zPad = [0,0];
zMask = [];

% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0;
params.wLam = 25;
params.gamma = 1;
params.maxIterReg = 800;
params.beta = 1.2;
params.maxIter = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;

P.params = params;

% coupled params
load([param_dir,'lambda_select_noise8'])
Pc.lambda_values = param_select;
Pc.initialization = 'simultaneous';
Pc.wLam = 25;
Pc.gamma = 1;
Pc.maxIterReg = 800;
Pc.num_outer_iters = 10;
Pc.baseFileName = 'fista_fit_%i_%i.mat';
Pc.num_ims = num_ims;
Pc.prefix = 'polar_vector';
Pc.dataset = [dataset];
% Gamma values
gamma_vals = [ 0.05,0.1, 0.2,0.01, 0.025, 0.035]; 
N = numel(gamma_vals);

% Select noise level
obj_gcv = zeros(N,num_ims);
err_gcv = zeros(N,num_ims);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
end
A0 = unshifted_basis_vector_stack_norm2(P);
%% Run grid search
for i = 6
    Pc.init_dir = [datadir,'noise_1D_coupled_simul_',dataset_num,'_init',num2str(i)];
    Pc.output_dirA = [datadir,'noise_1D_coupled_simul_',dataset_num,'_',num2str(i),'a'];
    Pc.output_dirB = [datadir,'noise_1D_coupled_simul_',dataset_num,'_',num2str(i),'b'];
    mkdir(Pc.init_dir)
    mkdir(Pc.output_dirA)
    mkdir(Pc.output_dirB)
    Pc.gamma = gamma_vals(i);
    runCoupledFISTA_1D(P,Pc)
end
% save([param_dir,'param_search_coupled_discrep_8.mat'],'err_gcv','obj_gcv','P','lambda_vals')

%% Select lambda values
N=6;
err_discrep = zeros(num_ims,N);
noise_eta = zeros(num_ims,N);
total_norm = zeros(N,1);
total_err = zeros(N,1);
for i = 1:N
    for image_num = 1:num_ims
        im_data = load([datadir,'noise_1D_coupled_',dataset_num,'_',num2str(i),'a\',sprintf(Pc.baseFileName,1,image_num)]);
        err_discrep(image_num,i) = im_data.err(end-1);
        b = im_data.polar_image;
        % Scale image by 2-norm
        bn = b/norm(b(:));

        %kernel = [0.0003 0.1329 0.9822 0.1329 0.0003];
        kernel = [0.0001 0.0003 0.0012 0.0042 0.0127 0.0329,...
                  0.0740 0.1434 0.2399 0.3466 0.4322 0.4652];
        kernel = [kernel,fliplr(kernel(1:end-1))];
        kernel = kernel./sum(kernel(:));

        bn_hat = conv(bn,kernel,'same');
        noise_eta(image_num,i) = norm(bn-bn_hat);
    end
    total_norm(i) = total_norm(i) + norm(b(:)).^2;
    total_err(i) = norm(err_discrep(:,i).^2);
    total_noise_est = 
end

gamma_index = find(discrep_crit == min(discrep_crit));


%% Compute true awmv
load(['D:\CHESS_data\simulated_two_spot_1D_noise2_12\synth_data.mat'])
truth_awmv_az = zeros(num_ims,1);

for i = 1:num_ims
    sample = synth_sample{i};
    az_awmv = sample.std_theta'*sample.amplitudes/sum(sample.amplitudes(:));
    truth_awmv_az(i) = az_awmv;
end

%% Load fit for selected parameters and plot
figure(222)
[ha1, pos1] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az = zeros(num_ims,1);
for image_num = 1:num_ims
    im_data = load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']));
    % Zero pad image
    b = im_data.polar_vector;
    % Scale image by 2-norm
    bn = b/norm(b(:));
    
    load(fullfile([datadir,'noise_1D_coupled_1a'],...
         sprintf(Pc.baseFileName,1,image_num)))
    
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    var_signal = squeeze(sum(x_hat,1));
    var_sum = sum(var_signal(:));
    awmv_az(image_num) = sum(sqrt(P.var_theta(:)).*var_signal(:))/var_sum;
    
    % Plot
    axes(ha1(image_num))
    hold on
    
    plot(bn)
    plot(fit)
    
    % load independent
    load(fullfile([datadir,'noise_1D_indep_param_8'],...
     sprintf(Pc.baseFileName,1,image_num)))
    
    var_signal = squeeze(sum(x_hat,1));
    var_sum = sum(var_signal(:));
    awmv_az2(image_num) = sum(sqrt(P.var_theta(:)).*var_signal(:))/var_sum;
    

end

figure(1)
hold on
plot(awmv_az)
plot(awmv_az2)
plot(truth_awmv_az)

% Also load up 

legend('coupled','indep','truth')
