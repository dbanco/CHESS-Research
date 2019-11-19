clear all
close all
P.set = 1;
dataset = 'D:\CHESS_data\simulated_two_spot_1D_noise2';
num_ims = 10;

%% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
load(fullfile([dataset,'_1'],[prefix,'_1.mat']));
P.num_theta= size(polar_vector,1);
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

%% fista params
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.t_k = 1;
params.lambda = 0.08;
params.wLam = 25;
params.beta = 1.2;
params.maxIter = 800;
params.maxIterReg = 800;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;
    


baseFileName = 'fista_fit_%i_%i.mat';

output_dir = 'D:\CHESS_data\simulated_two_spot_1D_noise2_independent';

lambda_vals = logspace(-2,0,10); 
N = numel(lambda_vals);
% noise level (1-12)
k = 8;
dir_num = sprintf('_%i',k);
obj_lambdas = zeros(N,num_ims);
err_lambdas = zeros(N,num_ims);
% lambda_vals = logspace(-3,0,500); 

for ii = 1:N
    P.params.lambda = lambda_vals(ii);
    % iterate over each image
    for image_num = 1:num_ims
        im_data = load(fullfile([dataset,dir_num],[prefix,'_',num2str(image_num),'.mat']));
        %% Zero pad image
        b = im_data.polar_vector;
        % Scale image by 2-norm
        b = b/norm(b(:));
    
        % Construct dictionary
        switch P.basis
            case 'norm2'
                A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
        end

        x_init = zeros(size(A0ft_stack));
        for i = 1:P.num_var_t
            x_init(:,i) = b/P.num_var_t;
        end

        [x_hat,err,obj,~,~,~] = FISTA_Circulant_1D(A0ft_stack,b,x_init,P.params);
        obj_lambdas(ii,image_num) = obj(end-1);
        err_lambdas(ii,image_num) = err(end-1);
        l1_term(ii,image_num) = norm(x_hat(:),1);
        fit = Ax_ft_1D(A0ft_stack,x_hat);
        res_term(ii,image_num) = 0.5*norm(b-fit)^2;
    end
end

save('param_search5.mat','err_lambdas','obj_lambdas','l1_term','res_term','P','lambda_vals')
%% Plot parameter search curves
load('param_search4.mat')
% lambda_vals = logspace(-3,0,500)
close all

figure(1)
plot(lambda_vals,err_lambdas)

figure(2)
plot(lambda_vals,obj_lambdas)

figure(3)
plot(res_term,l1_term)
xlabel('Residual')
ylabel('L-1')
title('L-curve')

figure(4)
for i = 1:num_ims
    subplot(5,2,i)
    hold on
    plot(log(res_term(:,i)),log(l1_term(:,i)),'-o')
    xlabel('log-Residual')
    ylabel('log-L-1')
end

image_num = 10;
xx = res_term(:,image_num)

min(find(xx>0.07607))
figure(5)
loglog(res_term(:,image_num),l1_term(:,image_num),'o-')
xlabel('Residual')
ylabel('L-1')

% figure(6)
% plot(res_term(:,image_num),l1_term(:,image_num),'o-')
% xlabel('Residual')
% ylabel('L-1')


% X = [res_term(:,noise_val),l1_term(:,noise_val)];
% [L,R,k] = curvature(X);
% figure(6)
% plot(k2)