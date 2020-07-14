clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;
% Parent directory
top_dir = 'D:\MMPAD_data';
%     top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring1_zero';

% Output dirs
output_name = '_indep_ISM1';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);

baseFileName = 'indep_fit_%i_%i.mat';

% Load most parameters by loading single output
load(fullfile(output_dir,sprintf(baseFileName,1,1)))
N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;
M = numel(P.lambda_values);

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);


%% Select lambda values
disp('Selecting lambda values')

% VDF figure
vdf_fig = figure(1);
[ha1, pos1] = tight_subplot(3,ceil(M/3),[0.1 0.03],[.02 .08],[.02 .02]); 
set(vdf_fig, 'Position',  [100, 100, 1000, 400])

vdf_time = zeros(M,T,K);
err_select = zeros(M,T);
l0_select = zeros(M,T);
l1_select = zeros(M,T);
x_indep = cell(T,1);
tv_time = zeros(M,T-1);
im_ind = 1;
for i = 1:M
    fprintf('%i of %i \n',i,M)
    for j = 1:T
        b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
        e_data = load(fullfile(output_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,e_data.x_hat),129:133);
        b = P.dataScale*sum(b_data.polar_image,1);
        err_select(i,j) = sum(( fit(:) - b(:) ).^2);
        l0_select(i,j) = sum(e_data.x_hat(:) > 1e-4*sum(e_data.x_hat(:)));
        l1_select(i,j) = sum(e_data.x_hat(:));
        az_signal = squeeze(sum(e_data.x_hat,1));
        var_sum = sum(az_signal(:));
        vdf_time(i,j,:) = az_signal/var_sum;
    end
    
%     axes(ha1(im_ind))
%     imagesc(squeeze(vdf_time(i,:,:)))
%     shading interp
%     caxis([0 0.6])
%     colormap(jet)
    
    title(['\lambda = ',sprintf('%1.1d',P.lambda_values(i))])
%     ylabel('t')
%     xlabel('\sigma')
    im_ind = im_ind + 1;
end
err_select(err_select > 10^10) = 0;
l0_select(l0_select > 10^10) = 0;
l1_select(l1_select > 10^10) = 0;

%% Criterion separate params
noise_eta = 0.2;
discrep_crit = abs(err_select'-noise_eta);

[lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
param_select = P.lambda_values(lambda_indices);
lambda_values_separate = param_select;

% Criterion single param
discrep_crit = abs(mean(err_select,2)-noise_eta);
lambda_index = find(discrep_crit == min(discrep_crit));
param_select_single = P.lambda_values(lambda_index);
lambda_values_single = ones(T,1)*param_select_single;

figure(2)
plot(param_select,'o-')
title('Parameters selected')


%% Plot
lambda_vals = P.lambda_values;

figure(3)
semilogx(lambda_vals,mean(l1_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_1 term')

% Plot
figure(4)
semilogx(lambda_vals,mean(err_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('error')

% Plot
figure(5)
semilogx(lambda_vals,mean(l0_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_0')


% Analyze L-curve
total_err = sum(err_select(1:M,:),2);
total_l1 = sum(l1_select(1:M,:),2);
total_err = total_err./max(total_err);
total_l1 = total_l1./max(total_l1);

% Find kink in L-cureve method #1
% slopes = (total_l1(2:end) - total_l1(1:end-1))./...
%      (total_err(2:end) - total_err(1:end-1));
% slope_select = find((slopes)>0);
% select_ind = slope_select(1);

% Find kink in L-cureve method #2
sq_origin_dist = total_l1.^2 + total_err.^2;
select_ind = find(sq_origin_dist == min(sq_origin_dist));

figure(6)
plot(total_l1,total_err,'o-')
% plot(mean(l1_select,2),mean(err_select,2),'o-')
hold on
plot(total_l1(select_ind),total_err(select_ind),'s','Markersize',15)
xlabel('time-average l1-norm')
ylabel('time-average error')
title('L-curve')

%% L curve parameter selection
select_indices = zeros(T,1);
for t = 1:T
    err_t = err_select(1:M,t);
    l1_t = l1_select(1:M,t);
    err_t = err_t;%/max(err_t);
    l1_t = l1_t;%/max(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist + (err_t == 0) )  );
end

for t = 1:T
    load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    b = P.dataScale*sum(polar_image,1);
    rel_err_t = err_select(1:M,t)/sum(b(:).^2);
    while rel_err_t(select_indices(t)) > 0.02
        if select_indices(t) > 1
            select_indices(t) = select_indices(t) - 1;
        else
            select_indices(t) = find(rel_err_t == min(rel_err_t));
            break
        end
    end
end

% figure(7)
% [ha7, ~] = tight_subplot(10,11,[.05 .05],[.01 .01],[.01 .01]); 
% tt = 1;
% for t = 1:5:T
%     err_t = err_select(1:M,t);
%     l1_t = l1_select(1:M,t);
%     err_t = err_t;%/max(err_t);
%     l1_t = l1_t;%/max(l1_t);
%     axes(ha7(tt))
%     plot(l1_t,err_t,'o-')
%     hold on
%     plot(l1_t(select_indices(t)),err_t(select_indices(t)),'s','MarkerSize',15)
%     hold on
%     tt = tt + 1;
% end

% xlabel('l1-norm')
% ylabel('error')
% title('L-curve')

%% Selected VDF;
select_vdf = zeros(T,K);
for j = 1:T
    select_vdf(j,:) = vdf_time(select_indices(j),j,:);
end
figure(8)
imagesc(select_vdf')
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

%% Plot fits for L-curve selected parameter
fits_fig = figure(9);
[ha2, ~] = tight_subplot(10,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
for t = 201:300
    load(fullfile(output_dir,sprintf(baseFileName,select_indices(t),t)))
    load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    
    polar_vector = sum(polar_image,1);
    fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_hat),129:133);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = P.dataScale*zeroPad(polar_vector,P.params.zeroPad);
    
    final_thresh = 1e-3*sum(x_hat(:));
    x_hat(x_hat<final_thresh) = 0;
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    rel_err = sum((fit(:)-b(:)).^2)./sum(b(:).^2);
    legend(sprintf('%i \n %0.2f',sum(x_hat(:)>0),rel_err),'location','northeast')
    im_ind = im_ind + 1;
end

%% Show coefficients of selected parameters
fits_fig = figure(10);
[ha2, ~] = tight_subplot(2,T/2,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
for t = 1:T
    load(fullfile(output_dir,sprintf(baseFileName,select_indices(t),t)))
    axes(ha2(t))
    imagesc(x_hat)
end

% %% Plot fits single selected paramter
% fits_fig = figure(223);
% [ha3, ~] = tight_subplot(2,T/2,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az_vdfs = zeros(T,1);
% im_ind = 1;
% for image_num = 1:T
% 
%     load(fullfile(output_dir,sprintf(baseFileName,lambda_index,image_num)))
%     load(fullfile(dataset,[P.prefix,'_',num2str(image_num),'.mat']))
%     
%     polar_vector = sum(polar_image,1);
%     fit = Ax_ft_1D(A0ft_stack,x_hat);
%     az_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(az_signal(:));
%     awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     b = zeroPad(polar_vector,P.params.zeroPad);
%     
%     final_thresh = 1e-4*sum(x_hat(:));
%     x_hat(x_hat<final_thresh) = 0;
%     % Plot
%     axes(ha3(im_ind))
%     hold on
%     plot(b)
%     plot(fit)
%     legend(sprintf('%i',sum(x_hat(:)>0)),'location','northeast')
%     im_ind = im_ind + 1;
% end
% 
% 
% %% Plot fits separate paramaters
% fits_fig = figure(224);
% [ha2, pos2] = tight_subplot(2,T/2,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az_vdfs = zeros(T,1);
% im_ind = 1;
% for image_num = 1:T
% 
%     load(fullfile(output_dir,sprintf(baseFileName,19,image_num)))
%     load(fullfile(dataset,[P.prefix,'_',num2str(image_num),'.mat']))
%     
%     polar_vector = sum(polar_image,1);
%     fit = Ax_ft_1D(A0ft_stack,x_hat);
%     az_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(az_signal(:));
%     awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     b = zeroPad(polar_vector,P.params.zeroPad);
%     
%     final_thresh = 1e-4*sum(x_hat(:));
%     x_hat(x_hat<final_thresh) = 0;
%     % Plot
%     axes(ha2(im_ind))
%     hold on
%     plot(b)
%     plot(fit)
%     legend(sprintf('%i',sum(x_hat(:)>0)),'location','northeast')
%     im_ind = im_ind + 1;
% end