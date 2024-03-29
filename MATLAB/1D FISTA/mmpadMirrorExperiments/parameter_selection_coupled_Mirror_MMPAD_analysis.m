clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;
% Parent directory
top_dir = 'D:\MMPAD_data_nr1';
% top_dir = 'E:\PureTiRD_full';
%     top_dir = '/cluster/shared/dbanco02';

for ring_num = 2
% Input dirs
dset_name = sprintf('ring%i_zero',ring_num);

% Output dirs
output_name = '_coupled_CG_TVphi_Mirror5';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);

baseFileName = 'coupled_fit_%i.mat';

% Load parameters by loading single output
load(fullfile(output_dir,sprintf(baseFileName,1)))
[~,indepDir] = fileparts(P.indepDir);
indepName = P.indepName;

% Real lambda values
lambda2_values = P.lambda2_values;
% lambda_values = [logspace(-7,-5.1,10) logspace(-5,1,30)];
% P.lambda_values = lambda_values;

[N,K,T] = size(X_hat);
M = numel(lambda2_values);

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);


%% Select lambda values
disp('Selecting lambda values')

err_select = zeros(M,T);
l0_select = zeros(M,T);
l1_select = zeros(M,T);
tv_penalty = zeros(M,1);
x_indep = cell(T,1);
vdfs = zeros(K,M,T);
awmv_az = zeros(M,T);

B = zeros(N,T);
% Load data
for j = 1:T
    b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
    b = P.dataScale*sum(b_data.polar_image,1);
    
    % Mirror data
    nn = numel(b);
    pad1 = floor(nn/2);
    pad2 = ceil(nn/2);
    N = nn + pad1 + pad2;
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+nn)) = b;
    b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flip(b(1:pad1));
    B(:,j) = flip(b_mirror);
end

tv_time = zeros(M,T-1);
im_ind = 1;
for i = 1:M
    fprintf('%i of %i \n',i,M)
    e_data = load(fullfile(output_dir,sprintf(baseFileName,i)),'P','X_hat');
    
    tv_penalty(i) = sum(abs(DiffPhiX_1D(e_data.X_hat)),'all');
    
    for j = 1:T
        b = B(:,j);
        x = squeeze(e_data.X_hat(:,:,j));
        fit = Ax_ft_1D(A0ft_stack,x);   
        err_select(i,j) = sum(( fit(:) - b(:) ).^2)/norm(b)^2;
        l0_select(i,j) = sum(x(:) > 1e-4*sum(x(:)));
        l1_select(i,j) = sum(x(:));
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        vdfs(:,i,j) = az_signal./var_sum;
        awmv_az(i,j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end

    im_ind = im_ind + 1;
end
err_select(err_select > 10^10) = 0;
l0_select(l0_select > 10^10) = 0;
l1_select(l1_select > 10^10) = 0;

% Load statistics for independently fit data
awmv_az_init = zeros(T,1);
err_indep = zeros(T,1);
l0_indep = zeros(T,1);
l1_indep = zeros(T,1);
vdfs_indep = zeros(P.num_var_t,T);
X_indep = zeros(N,K,T);

for j = 1:T
    i_data = load(fullfile(top_dir,indepDir,sprintf(indepName,e_data.P.params.lambda1_indices(j),j)),'err','x_hat');
    load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']) )
    x = i_data.x_hat;
    X_indep(:,:,j) = x;
    fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x),P.params.zeroMask);
    b = B(:,j);
    err_indep(j) = sum((fit(:)-b(:)).^2)/sum(b(:).^2);
    l0_indep(j) = sum(x(:)>0);
    l1_indep(j) = sum(x(:));
    az_signal = squeeze(sum(x,1));
    var_sum = sum(az_signal(:));
    vdfs_indep(:,j) = az_signal./var_sum;
    awmv_az_init(j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end

tv_indep = sum(abs(DiffPhiX_1D(X_indep)),'all');

%% Plot awmv
figure_dir = ['C:\Users\dpqb1\OneDrive\Desktop\',output_subdir,'_figures\'];
mkdir(figure_dir)

% Plot AWMV
awmv_fig = figure(1);
legend_str = {};
% legend_str{1} = 'truth';
legend_str{1} = 'indep';
kk = 2;
hold on
plot(awmv_az_init,'LineWidth',1.5)
for k = 25
    hold on
    plot(awmv_az(k,:),'LineWidth',1.5)
    legend_str{kk} = sprintf('%0.3f',lambda2_values(k));
    kk = kk + 1;
end
% ylim([0 40])
ylabel('AWMV_\eta','FontSize',20)
xlabel('t','FontSize',20)
legend(legend_str,'location','best')
% save([dset_name,'_mirror_coupled_awmv.mat'],'awmv_az','awmv_az_init','lambda2_values')

%% Plot fits for L-curve selected parameter
load(fullfile(output_dir,sprintf(baseFileName,25)))
fits_fig = figure(99);
[ha2, ~] = tight_subplot(5,5,[.005 .005],[.01 .01],[.01 .01]); 
im_ind = 1;
for t = 166:190
    x_hat = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    b = B(:,t);

    % Plot
    axes(ha2(im_ind))
    hold on
    plot((b))
    plot(fit)
    rel_err = sum((fit(:)-b(:)).^2)/norm(b(:))^2;
    legend(sprintf('%i',t),'location','northeast')
    im_ind = im_ind + 1;
end

%% Plot fits to show in figure
load(fullfile(output_dir,sprintf(baseFileName,25)))
fits_fig = figure(99);
[ha2, ~] = tight_subplot(1,5,[.005 .005],[.01 .01],[.01 .01]); 
im_ind = 1;
for t = [1,40,80,130,240]
    x_hat = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    b = B(:,t);

    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b((pad1+1):(pad1+nn)),'Linewidth',2)
    plot(fit((pad1+1):(pad1+nn)),'Linewidth',1)
    rel_err = sum((fit(:)-b(:)).^2)/norm(b(:))^2;
    im_ind = im_ind + 1;
    set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'XColor', 'none','YColor','none')
end
%     legend('data','fit','location','best')
[1,40,80,130,240]/4
%% Plot indep fits to show in figure
load(fullfile(output_dir,sprintf(baseFileName,25)))
fits_fig = figure(96);
[ha2, ~] = tight_subplot(1,5,[.005 .005],[.01 .01],[.01 .01]); 
im_ind = 1;
for t = [1,40,80,130,240]
    x_hat = X_indep(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    b = B(:,t);

    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b((pad1+1):(pad1+nn)),'Linewidth',2)
    plot(fit((pad1+1):(pad1+nn)),'Linewidth',1)
    rel_err = sum((fit(:)-b(:)).^2)/norm(b(:))^2;
    im_ind = im_ind + 1;
    set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'XColor', 'none','YColor','none')
end
%     legend('data','fit','location','best')
[1,40,80,130,240]/4
%% Plot pointwiseawmv
load(fullfile(output_dir,sprintf(baseFileName,25)))
pw_fig = figure(97);
[hapw, ~] = tight_subplot(5,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
pwawmv = zeros(N,1);
for t = 166:190
    x_hat = X_hat(:,:,t);
    for i = 1:N
        pwawmv(i) = sum( sqrt(P.var_theta(:)).*x_hat(i,:)' );
    end

    % Plot
    axes(hapw(im_ind))
    hold on
    plot(pwawmv)
    legend(sprintf('%i',t),'location','northeast')
    im_ind = im_ind + 1;
end


end

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
semilogx(lambda_vals(1:M),mean(l1_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('l_1 term')

% Plot
figure(4)
semilogx(lambda_vals(1:M),mean(err_select,2),'o-')
hold on
xlabel('\lambda')
ylabel('error')

% Plot
figure(5)
semilogx(lambda_vals(1:M),mean(l0_select,2),'o-')
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

figure(666)
iii=55;
loglog(l1_select(1:M,iii),err_select(1:M,iii),'o-')
hold on
loglog(l1_select(select_ind,iii),err_select(select_ind,iii),'s','Markersize',15)
xlabel('l1-norm')
ylabel('error')
title('L-curve')

figure(667)
plot(l1_select(1:M,iii),err_select(1:M,iii),'o-')
hold on
plot(l1_select(select_ind,iii),err_select(select_ind,iii),'s','Markersize',15)
xlabel('l1-norm')
ylabel('error')
title('L-curve')

%% L curve parameter selection
select_indices = zeros(T,1)+1;
for t = 1:T
    err_t = err_select(1:M,t);
    l1_t = l1_select(1:M,t);
    err_t = err_t;%/max(err_t);
    l1_t = l1_t/max(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist + (err_t == 0) )  );
    if t == 128
        figure(1111)
        hold on
        plot(l1_t,err_t,'o-')
        plot(l1_t(select_indices(t)),err_t(select_indices(t)),'s')
        xlabel('l1-norm')
        ylabel('error')
    end
    
    
end

for t = 1:T
    b = B(:,t);
    rel_err_t = err_select(1:M,t);
    while rel_err_t(select_indices(t)) > 0.02
        if select_indices(t) > 1
            select_indices(t) = select_indices(t) - 1;
        else
            select_indices(t) = find(rel_err_t == min(rel_err_t));
            break
        end
    end
end

figure(543)
plot(select_indices)

%% Plot AWMV
load(fullfile(output_dir,sprintf(baseFileName,1)))
for t = 1:T
%     load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    x_hat = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end
figure(12)
plot(awmv_az_vdfs)

%% Selected VDF;
select_vdf = zeros(T,K);
for j = 1:T
    select_vdf(j,:) = vdfs(select_indices(j),j,:);
end
figure(8)
imagesc(select_vdf')
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')

%% Plot fits for L-curve selected parameter
load(fullfile(output_dir,sprintf(baseFileName,1)))
fits_fig = figure(9);
[ha2, ~] = tight_subplot(5,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
for t = 1:5:200 
    x_hat = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = B(:,t);

    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    rel_err = sum((fit(:)-b(:)).^2)/norm(b(:))^2;
    legend(sprintf('%i \n %0.2f',sum(x_hat(:)>0),rel_err),'location','northeast')
    im_ind = im_ind + 1;
end


%}