clear all
close all

num_theta = 180;
noise_std = 0:0.03:0.30;
n_eta_levels = sqrt(num_theta.*noise_std.^2);
% n_eta_levels = noise_std;
% n_eta_levels = linspace(0.02,0.35,numel(noise_std));

n_level = 3;

% Parameter selection
disp('Setup params')
P.set = 1;

dset_name = 'gnoise4_nonorm';
num_ims = 20;

% datadir = '/cluster/shared/dbanco02/';
% dataset = ['/cluster/home/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'/'];
% indep_dir = ['/cluster/shared/dbanco02/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep_approx2/'];
% init_dir = [datadir,'gnoise4_subdir/simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init_indep_approx2'];

datadir = 'D:\CHESS_data\';
dataset = ['D:\CHESS_data\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'\'];
indep_dir = ['D:\CHESS_data\ADMM_CG_indep1\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_indep\'];
init_dir = ['D:\CHESS_data\ADMM_CG_indep1\simulated_two_spot_1D_',dset_name,'_',num2str(n_level),'_simul_init_indep'];

% Universal Parameters
% Ring sampling parameters
prefix = 'polar_vector';
baseFileName = 'fista_fit_%i.mat';

% Lambda values
lambda_vals = logspace(-4,1,30); 
N = numel(lambda_vals);

% Load most parameters by loading single output
load([indep_dir,sprintf(baseFileName,1)])

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
    x_data = load(fullfile(indep_dir,sprintf(baseFileName,i)));

    for j = 1:num_ims
        x = x_data.X_hat(:,:,j);
        fit = Ax_ft_1D(A0ft_stack,x);
        bb = x_data.polar_image;
        err_select(i,j) = sum( (fit(:)-bb(:)).^2);
        l0_select(i,j) = sum(x(:) > 1e-4*sum(x(:)));
        l1_select(i,j) = sum(x(:));
        az_signal = squeeze(sum(x,1));
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
figure(2)
semilogx(lambda_vals,mean(l1_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_1 term')

figure(3)
semilogx(lambda_vals,mean(err_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('error')

figure(4)
semilogx(lambda_vals,mean(l0_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_0')

figure(5)
total_err = sum(err_select(1:20,:),2);
total_l1 = sum(l1_select(1:20,:),2);
plot(total_err,total_l1,'o-')
xlabel('Error')
ylabel('l-1 norm')
title('L-curve')

% select parameter from L-curve
slopes = (total_l1(2:end) - total_l1(1:end-1))./...
         (total_err(2:end) - total_err(1:end-1));
slope_select = find(abs(slopes)<1);
slope_select = slope_select(1);
selectx = total_err(slope_select);
selecty = total_l1(slope_select);
hold on
plot(selectx,selecty,'s','Markersize',14)

% figure(545)
% plot(total_err(1:end-1),slopes,'o-')
% title('slope')

figure(544)
plot((curvature),'o-')
title('curvature')
% figure(523)
% loglog(coord(:,1),coord(:,2),'o-')

%% Plot fits separate selected paramters
fits_fig = figure(222);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
az_vdf = zeros(num_ims,P.num_var_t);
im_ind = 1;
trial_k = 1;
for image_num = 1:20
    load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']) )
    
    x_data = load(fullfile(indep_dir,sprintf(baseFileName,lambda_indices(image_num))));
    x = x_data.X_hat(:,:,image_num);
    
    az_signal = squeeze(sum(x,1));
    var_sum = sum(az_signal(:));
    az_vdf(image_num,:) = az_signal./var_sum;
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x(:));
    x(x<final_thresh) = 0;
    fit = Ax_ft_1D(A0ft_stack,x);
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end

figure(2221)
subplot(2,1,1)
imagesc(az_vdf)
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

subplot(2,1,2)
plot(awmv_az_vdfs,'o-')
ylabel('AWMV')
xlabel('time')
%% Plot fits single selected paramter
fits_fig = figure(223);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
im_ind = 1;
trial_k = 1;
for image_num = 1:20
    load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']) )
    
    x_data = load(fullfile(indep_dir,sprintf(baseFileName,lambda_index)));
    x = x_data.X_hat(:,:,image_num);
    
    az_signal = squeeze(sum(x,1));
    var_sum = sum(az_signal(:));
    az_vdf(image_num,:) = az_signal./var_sum;
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x(:));
    x(x<final_thresh) = 0;
    fit = Ax_ft_1D(A0ft_stack,x);
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end

figure(2231)
subplot(2,1,1)
imagesc(az_vdf)
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

subplot(2,1,2)
plot(awmv_az_vdfs,'o-')
ylabel('AWMV')
xlabel('time')

%% Plot fits with a L curve paramater
fits_fig = figure(224);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
im_ind = 1;
trial_k = 1;
for image_num = 1:20
    load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']) )
    
    x_data = load(fullfile(indep_dir,sprintf(baseFileName,L_select_ind)));
    x = x_data.X_hat(:,:,image_num);
    
    az_signal = squeeze(sum(x,1));
    var_sum = sum(az_signal(:));
    az_vdf(image_num,:) = az_signal./var_sum;
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x(:));
    x(x<final_thresh) = 0;
    fit = Ax_ft_1D(A0ft_stack,x);
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end

figure(2241)
subplot(2,1,1)
imagesc(az_vdf)
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

subplot(2,1,2)
plot(awmv_az_vdfs,'o-')
ylabel('AWMV')
xlabel('time')

%% Plot fits with a specified paramater
fits_fig = figure(224);
[ha2, pos2] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(num_ims,1);
im_ind = 1;
trial_k = 1;
param_index = curve_select;
for image_num = 1:20
    load(fullfile(dataset,[prefix,'_',num2str(image_num),'.mat']) )
    
    x_data = load(fullfile(indep_dir,sprintf(baseFileName,param_index)));
    x = x_data.X_hat(:,:,image_num);
    
    az_signal = squeeze(sum(x,1));
    var_sum = sum(az_signal(:));
    az_vdf(image_num,:) = az_signal./var_sum;
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-4*sum(x(:));
    x(x<final_thresh) = 0;
    fit = Ax_ft_1D(A0ft_stack,x);
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x(:)>0)),'location','northeast')
    im_ind = im_ind + 1;
end

figure(2241)
subplot(2,1,1)
imagesc(az_vdf)
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

subplot(2,1,2)
plot(awmv_az_vdfs,'o-')
ylabel('AWMV')
xlabel('time')
