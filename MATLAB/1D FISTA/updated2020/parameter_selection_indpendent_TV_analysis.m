clear all
close all

num_theta = 180;
noise_std = 0:0.03:0.30;
n_eta_levels = sqrt(num_theta.*noise_std.^2);
% n_eta_levels = noise_std;
% n_eta_levels = linspace(0.02,0.35,numel(noise_std));

n_level = 3;

% Parameter selection
disp('Setup parms')
P.set = 1;

dset_name = 'gnoise4_nonorm';
num_ims = 20;

% datadir = '/cluster/shared/dbanco02/';
% dataset = ['/cluster/home/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'/'];
% indep_dir = ['/cluster/shared/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep_approx2/'];
% init_dir = [datadir,'gnoise4_subdir/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init_indep_approx2'];

datadir = 'D:\CHESS_data\';
dataset = ['D:\CHESS_data\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'\'];
indep_dir = ['D:\CHESS_data\ADMM_CG_indep1\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep1\'];
init_dir = ['D:\CHESS_data\ADMM_CG_indep1\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init_indep1'];

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
baseFileName = 'fista_fit_%i_%i.mat';

% Lambda values
lambda_vals = logspace(-4,1,30); 
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
x_indep = cell(num_ims,1);
tv_time = zeros(N,num_ims-1);
im_ind = 1;
for i = 1:N
    fprintf('%i of %i \n',i,N)
    for j = 1:num_ims
        e_data = load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        err_select(i,j) = e_data.err(end);
        l0_select(i,j) = sum(e_data.x_hat(:) > 1e-4*sum(e_data.x_hat(:)));
        l1_select(i,j) = sum(e_data.x_hat(:));
        az_signal = squeeze(sum(e_data.x_hat,1));
        var_sum = sum(az_signal(:));
        vdf_time(i,j,:) = az_signal/var_sum;
        
%         x_indep{j} = e_data.x_hat;
    end
%     for j = 1:num_ims-1
%         tv_dist = sum( sqrt( ( x_indep{j}(:) - x_indep{j+1}(:) ).^2 + tvBeta^2) );
%         tv_time(i,j) = tv_dist;
%     end
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

% Criterion separate params
noise_eta = n_eta_levels(n_level);
discrep_crit = abs(err_select'-noise_eta);

[lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
param_select = lambda_vals(lambda_indices);
lambda_values_separate = param_select;

% Criterion single param
discrep_crit = abs(mean(err_select,2)-noise_eta);
lambda_index = find(discrep_crit == min(discrep_crit));
param_select_single = lambda_vals(lambda_index);
lambda_values_single = ones(num_ims,1)*param_select_single;

figure(112)
plot(param_select,'o-')
title('Parameters selected')

%% Selected VDF
select_vdf = zeros(num_ims,size(A0,2));
for i = 1:N
    for j = 1:num_ims
        select_vdf(j,:) = vdf_time(lambda_indices(j),j,:);
    end
end
figure(7)
imagesc(select_vdf)
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

%% Plot
for image_num = 3;
figure(2)
semilogx(lambda_vals,mean(l1_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_1 term')
end

% Plot
for image_num = 3;
figure(3)
semilogx(lambda_vals,mean(err_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('error')
end

% % Plot
% for image_num = 3;
% figure(6)
% plot(mean(l1_approx,2),mean(err_select,2),'o-')
% hold on
% xlabel('time-average l1-norm')
% ylabel('time-average error')
% title('L-curve')
% end


% Plot
for image_num = 3;
figure(4)
semilogx(lambda_vals,mean(l0_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_0')
end

%% Plot
% for image_num = 3;
% figure(2)
% semilogx(lambda_vals,l1_approx(:,image_num),'o-')
% hold on
% xlabel('\lambda')
% ylabel('l_1 term')
% end
% 
% % Plot
% for image_num = 3;
% figure(3)
% semilogx(lambda_vals,err_select(:,image_num),'o-')
% hold on
% xlabel('\lambda')
% ylabel('error')
% end
% 
% % Plot
% for image_num = 3;
% figure(4)
% semilogx(lambda_vals,l0_select(:,image_num),'o-')
% hold on
% xlabel('\lambda')
% ylabel('l_0')
% end
%% Plot fits separate selected paramters
fits_fig = figure(222);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
im_ind = 1;
trial_k = 1;
for image_num = 1:20

    load(fullfile(indep_dir,sprintf(baseFileName,lambda_indices(image_num),image_num)))

    polar_vector = polar_image;
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x_hat(:));
    x_hat(x_hat<final_thresh) = 0;
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x_hat(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end

%% Plot fits single selected paramter
fits_fig = figure(223);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
im_ind = 1;
trial_k = 1;
for image_num = 1:20

    load(fullfile(indep_dir,sprintf(baseFileName,lambda_index,image_num)))

    polar_vector = polar_image;
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x_hat(:));
    x_hat(x_hat<final_thresh) = 0;
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x_hat(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end


%% Plot fits separate paramaters
fits_fig = figure(224);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
im_ind = 1;
trial_k = 1;
for image_num = 1:20

    load(fullfile(indep_dir,sprintf(baseFileName,21,image_num)))

    polar_vector = polar_image;
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x_hat(:));
    x_hat(x_hat<final_thresh) = 0;
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x_hat(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end

%% Plot AWMV
% figure(23)
% plot(awmv_az_vdfs)

%% View convergence
for i = 1:N
    fprintf('%i of %i \n',i,N)
    for j = 1:num_ims
        load(fullfile(indep_dir,sprintf(baseFileName,i,j)),'obj');
        figure(99)
        plot(obj)
        pause
    end
end