%% Problem parameters
% Data I/O directories
data_dir = '/cluster/home/dbanco02/al7075_311_polar/';
results_dir = 'D:\CHESS_results\test\';

% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 2048;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,pi/64,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  2,       P.num_var_r).^2;

% Generate unshifted basis function matrices
P.betap = P.dtheta*P.drad;
P.weight = 1;
P.alphap = 10;


A0ft_stack = unshifted_basis_matrix_ft_stack(P);
A0_stack = unshifted_basis_matrix_stack(P);

% FISTA parameters
params.stoppingCriterion = 2;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 10;
params.gamma = 0.5;
params.maxCycles = 10;
params.isNonnegative = 1;
P.params = params;

P.load_step = 0;
% Map img num for row,col coordinate on sample
N = 41;
M = 5;
row = N;
col = 1;
img_array = zeros(N,M);
for img = 0:204
    img_array(row,col) = img;
    col = col + 1;
    if col > M
        row = row - 1;
        col = 1;
    end
end

% Reduce to subset of sample (41x5 -> 3x5)
P.img_array = img_array;
[N,M] = size(P.img_array);

% Reduce size of each image to be 21x50
for i = 1:N
    for j = 1:M
        % Load polar_image
        load([data_dir,... 
        'polar_image_',...
        num2str(P.load_step),'_',...
        num2str(P.img_array(i,j)), '.mat']);

        % Reduce image 
        test_im = polar_image(:,:);
        b_array{i,j} = test_im;
    end
end
%% cyclic FISTA with backtracking
c = parcluster('local');
c.NumWorkers = 32;
parpool(c, c.NumWorkers);
[x_hat_array, error_array]  = cyclic_optimization(A0ft_stack,b_array,P.params);   
save('cyclic_fista_fit_load_0.mat','x_hat_array','params')

%% Show results
% var_signal = zeros(P.num_var_t,P.num_var_r,N,M);
% hvs_az = zeros(N,M);
% k = 1;
% for i = 1:N
%     for j = 1:M
%         figure(1)
%         subplot(3,5,k)
%         im_fit = Ax_ft_2D(A0ft_stack,x_hat_array{i,j});
%         image((im_fit))
%         colormap jet
%         
%         figure(2)
%         subplot(3,5,k)
%         image((b_array{i,j}))
%         colormap jet
%         
%         k = k + 1;
%         
%         var_signal(:,:,i,j) = squeeze(sum(sum(x_hat_array{i,j},1),2)); 
%         hvs_az(i,j) = squeeze(sum(sum(var_signal(6:end,:,i,j))))/sum(sum(var_signal(:,:,i,j)));
%     end
% end
% figure(3)
% imagesc(hvs_az)
% colormap jet
% colorbar()
