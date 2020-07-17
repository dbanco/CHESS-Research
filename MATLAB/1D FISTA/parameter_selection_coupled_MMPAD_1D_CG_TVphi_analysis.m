%% Parameter selection
clear all
close all

% Parameter selection
disp('Setup parms')
P.set = 1;
% Parent directory
top_dir = 'D:\MMPAD_data';

% Input dirs
dset_name = 'ring1_zero_subset';

% Indep dirs
indep_name = '_indep_ISM1';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
output_name = '_coupled_CG_TVphi1';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);

% File Parameters
prefix = 'mmpad_img';
baseFileName = 'coupled_fit_%i.mat';

% Universal Parameters
% Ring sampling parameters
load(fullfile(output_dir,sprintf(baseFileName,1)))
lambda2_vals = P.lambda2_values;
M = numel(lambda2_vals);

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
T = P.num_ims;
[N,K] = size(A0ft_stack);

for k = 1:M
    fprintf('%i of %i\n',k,M)
    x_data = load( fullfile(output_dir,sprintf(baseFileName,k)) );
    
    for j = 1:P.num_ims
        % Fit objective
        x = x_data.X_hat(:,:,j);
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x),P.params.zeroMask);

        load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']) )
        b = sum(polar_image,1);
        err_select(k,j) = sum( (b(:)-fit(:)).^2 ) ;

        % Vdf
        az_signal = squeeze(sum(x,1));
        var_sum = sum(az_signal(:));
        vdfs(:,k,j) = az_signal./var_sum;
        awmv_az(k,j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;

        % Sparsity objective
        l0_select(k,j) = sum(x(:)>0);
        l1_select(k,j) = sum(x(:));
    end
end

% Compute tv objective
tv_penalty = zeros(M,1);
for k = 1:M
    fprintf('TVx %i of %i\n',k,M)
    x_data = load(fullfile(output_dir,sprintf(P.baseFileName,k)));
    tv_penalty(k) = sum(abs(DiffPhiX_1D(x_data.X_hat)),'all');
end

% Load statistics for independently fit data
awmv_az_init = zeros(T,1);
err_indep = zeros(T,1);
l0_indep = zeros(T,1);
l1_indep = zeros(T,1);
vdfs_indep = zeros(P.num_var_t,T);
X_indep = zeros(N,K,T);
indepFileName = 'indep_fit_%i_%i.mat';
for j = 1:T
    e_data = load(fullfile(indep_dir,sprintf(indepFileName,x_data.P.params.lambda1_indices(j),j)),'err','x_hat');
    load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']) )
    x = e_data.x_hat;
    X_indep(:,:,j) = x;
    fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x),P.params.zeroMask);
    b = sum(polar_image,1);
    err_indep(j) = sum((fit(:)-b(:)).^2);
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
legend_str{1} = '0';
kk = 2;
hold on
plot(awmv_az_init,'LineWidth',1.5)
for k = 1:M
    hold on
    plot(awmv_az(k,:),'LineWidth',1.5)
    legend_str{kk} = sprintf('%0.01s',lambda2_vals(k));
    kk = kk + 1;
end
% ylim([0 40])
ylabel('AWMV_\eta','FontSize',20)
xlabel('t','FontSize',20)
legend(legend_str,'location','best','FontSize',16)
%% Save images
% Plot err
total_err = sum(err_select,2);
figure(2)
hold on
% plot(0,sum(err_indep),'o')
plot(lambda2_vals,total_err,'o-')
ylabel('Error')
xlabel('Coupling parameter')

% Plot l1 norm
total_l1 = sum(l1_select,2);
figure(3)
hold on
% plot(0,sum(l1_indep),'o')
plot(lambda2_vals,total_l1,'o-')
ylabel('l1-norm')
xlabel('Coupling parameter')

% Plot number nonzeros coefficients
mean_l0 = mean(l0_select,2);
figure(4)
% loglog(0,mean(l0_indep),'o')
hold on
plot(lambda2_vals,mean_l0,'o-')
ylabel('l0-norm')
xlabel('Coupling parameter')


%     % Plot AWMV
%     figure(5)
%     legend_str = {};
%     legend_str{1} = 'indep';
%     hold on
%     plot(awmv_az_init,'-','LineWidth',1.5)
%     kk = 2;
%     for k = [4,5,M]
%         hold on
%         plot(awmv_az(sort_i(k),:),'-','LineWidth',1.5)
%         legend_str{kk} = sprintf('%0.04f',sort_lambda2(k));
%         kk = kk + 1;
%     end
%     legend(legend_str,'Location','Best')

% % Plot L-curve 1
% figure(6)
% plot(obj2(sort_i),obj1(sort_i),'o-')
% xlabel('l1-norm')
% ylabel('Error')

% Plot L-curve 2

Lcurve_fig = figure(7);
plot(total_err,tv_penalty,'o-')
ylabel('TV')
xlabel('Error')

total_err = sum(err_select,2);
total_l1 = sum(tv_penalty,2);

% Find kink in L-cureve method #1
% slopes = (total_l1(2:end) - total_l1(1:end-1))./...
%      (total_err(2:end) - total_err(1:end-1));
% slope_select = find(abs(slopes)<1);
% select_ind = slope_select(1);
select_ind = 1
% Find kink in L-cureve method #2
%     miny = min(tv_penalty); minx = min(total_err); coord =
%     [total_err./minx,tv_penalty./miny]; origin_dist = coord-1;
%     [val,select_ind] = min(sum(origin_dist.^2,2));

selectx = total_err(select_ind);
selecty = tv_penalty(select_ind,:);
hold on
plot(selectx,selecty,'s','Markersize',14)
plot(sum(err_indep),tv_indep,'*','Markersize',14)

%% Plot paramter selected
select_fig = figure(11);
legend_str = {};
legend_str{1} = 'indep';
hold on
plot(awmv_az_init,'LineWidth',1.5)
kk = 2;
for k = [1]
    hold on
    plot(awmv_az(k,:),'LineWidth',1.5)
    legend_str{kk} = sprintf('%0.01s',lambda2_vals(k));
    kk = kk + 1;
end

ylabel('AWMV_\eta','FontSize',20)
xlabel('t','FontSize',20)
legend(legend_str,'location','best','FontSize',16)
%     saveas(select_fig,[figure_dir,'awmv_select_',dset_name,'_',dataset_num,'.png'])

%% Plot vdf surface
vdf_time_all_fig = figure(56);
[ha3, pos3] = tight_subplot(3,ceil(M/3)+1,[0.1 0.03],[.02 .08],[.02 .02]); 
im_ind = 1;
    % Plot surface
axes(ha3(im_ind))
imagesc(squeeze(vdfs_indep'))
shading interp
caxis([0 0.6])
colormap(jet)

title(['\lambda_2 = ','0'])
ylabel('t')
xlabel('\sigma')

im_ind = im_ind + 1;

for trial_k = 1:M
    fprintf('%i of %i\n',trial_k,M)
    x_data = load(fullfile(output_dir,sprintf(baseFileName,trial_k)));
    for image_num = 1:T
        x_hat = x_data.X_hat(:,:,image_num);
        az_signal = squeeze(sum(x_hat,1));
        var_sum = sum(az_signal(:));
        vdf_time(trial_k,image_num,:) = az_signal/var_sum;
    end

    % Plot surface
    axes(ha3(im_ind))
    imagesc(squeeze(vdf_time(trial_k,:,:)))
    shading interp
    caxis([0 0.6])
    colormap(jet)

    title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(trial_k))])
    ylabel('t')
    xlabel('\sigma')

    im_ind = im_ind + 1;
end
saveas(vdf_time_all_fig,[figure_dir,'vdf_time_all_',dset_name,'_',dataset_num,'.png'])

vdf_time_fig = figure(566);
% Plot surface
imagesc(squeeze(vdf_time(select_ind,:,:)))
shading interp
caxis([0 0.6])
colormap(jet)

title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(select_ind))])
ylabel('t')
xlabel('\sigma')
%     saveas(vdf_time_fig,[figure_dir,'vdf_time_select_',dset_name,'_',dataset_num,'.png'])

vdf_indep_fig = figure(567);
% Plot surface
    imagesc(squeeze(vdfs_indep'))
shading interp
caxis([0 0.6])
colormap(jet)

title(['\lambda_2 = ','0'])
ylabel('t')
xlabel('\sigma')



%% Plot fits
fits_fig = figure(222);
[ha2, pos2] = tight_subplot(10,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
x_data = load(fullfile(output_dir,sprintf(baseFileName,1)));
for image_num = 1:100
    x_hat = x_data.X_hat(:,:,image_num);
    load(fullfile(dataset,[P.prefix,'_',num2str(image_num),'.mat']) )
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    polar_vector = sum(polar_image,1);
    b = P.dataScale*zeroPad(polar_vector,P.params.zeroPad);

    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',sum(x_hat(:)>1e-4)),'location','northeast')
    im_ind = im_ind + 1;
end
%     saveas(fits_fig,[figure_dir,'fits_',dset_name,'_',dataset_num,'.png'])


%% Plot vdfs
% figure(333)
% [ha3, pos3] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az_vdfs = zeros(num_ims,1);
% im_ind = 1;
% trial_k = 5;
% for image_num = 1:10
%     
%     load(fullfile([output_dir,num2str(sort_i(trial_k)),'a\'],sprintf(baseFileName,1,image_num)))
%     
%     polar_vector = squeeze(sum(polar_image,1));
%     fit = Ax_ft_1D(A0ft_stack,x_hat)*norm(polar_vector(:));
%     az_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(az_signal(:));
%     awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     % Plot
%     axes(ha3(im_ind))
%     plot(az_signal/var_sum)
%     legend(sprintf('%i',sum(x_hat(:)>1e-6)),'location','northeast')
%     im_ind = im_ind + 1;
% end

