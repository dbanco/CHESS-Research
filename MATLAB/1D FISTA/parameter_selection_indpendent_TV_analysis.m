clear all
close all

noise_std = 0:0.03:0.30;
n_eta_levels = 0.5*sqrt(180.*noise_std.^2) + 0.10;
% n_eta_levels = linspace(0.02,0.35,numel(noise_std));

n_level = 4;

% Parameter selection
disp('Setup parms')
P.set = 1;

dset_name = 'gnoise4_nonorm';
num_ims = 20;

datadir = '/cluster/shared/dbanco02/';
dataset = ['/cluster/home/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'/'];
indep_dir = ['/cluster/shared/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep_approx2/'];
init_dir = [datadir,'gnoise4_subdir/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init_indep_approx2'];


% datadir = 'E:\CHESS_data\';
% dataset = ['E:\CHESS_data\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'\'];
% indep_dir = ['E:\CHESS_data\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep_approx2\'];
% init_dir = [datadir,'\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init_indep_approx2'];

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-3,1,30); 
N = numel(lambda_vals);

% Load most parameters by loading single output
load([indep_dir,sprintf(baseFileName,1,1)])

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
A0 = unshifted_basis_vector_stack_norm2_zpad(P);

%% Select lambda values
disp('Selecting lambda values')

surf_fig = figure(111);
[ha1, pos1] = tight_subplot(3,ceil(N/3),[0.1 0.03],[.02 .08],[.02 .02]); 
set(surf_fig, 'Position',  [100, 100, 1000, 400])

vdf_time = zeros(N,num_ims,size(A0,2));
err_select = zeros(N,num_ims);
l0_select = zeros(N,num_ims);
l1_select = zeros(N,num_ims);
im_ind = 1;
for i = 1:N
    fprintf('%i of %i \n',i,N)
    for j = 1:num_ims
        e_data = load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        err_select(i,j) = e_data.err(end);
        l0_select(i,j) = sum(e_data.x_hat(:) > 0);
        l1_select(i,j) = sum(e_data.x_hat(:));
    
        az_signal = squeeze(sum(e_data.x_hat,1));
        var_sum = sum(az_signal(:));
        vdf_time(i,j,:) = az_signal/var_sum;
    end
    axes(ha1(im_ind))
    imagesc(squeeze(vdf_time(i,:,:)))
    shading interp
    caxis([0 0.6])
    colormap(jet)
    
    title(['\lambda = ',sprintf('%1.1d',lambda_vals(i))])
%     ylabel('t')
%     xlabel('\sigma')
    im_ind = im_ind + 1;
end
err_select(err_select > 10^10) = 0;
l0_select(l0_select > 10^10) = 0;
l1_select(l1_select > 10^10) = 0;

% Criterion 
noise_eta = n_eta_levels(n_level);
discrep_crit = abs(err_select'-noise_eta);

[lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
param_select = lambda_vals(lambda_indices);
Pc.lambda_values = param_select;

figure(112)
plot(param_select,'o-')

%% Plot fits
% 
% fits_fig = figure(222);
% [ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az_vdfs = zeros(num_ims,1);
% im_ind = 1;
% trial_k = 1;
% for image_num = 1:20
% 
%     load(fullfile(init_dir,sprintf(baseFileName,1,image_num)))
% 
%     polar_vector = polar_image;
%     fit = Ax_ft_1D(A0ft_stack,x_hat);
%     az_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(az_signal(:));
%     awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     b = zeroPad(polar_vector,P.params.zeroPad);
% 
%     % Plot
%     axes(ha2(im_ind))
%     hold on
%     plot(b)
%     plot(fit)
%     legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
%     im_ind = im_ind + 1;
% end