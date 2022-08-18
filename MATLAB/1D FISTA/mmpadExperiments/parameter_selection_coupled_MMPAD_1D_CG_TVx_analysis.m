% Parameter selection
disp('Setup params')
P.set = 1;

% Parent directory
top_dir = 'D:\MMPAD_data';
%     top_dir = '/cluster/shared/dbanco02';

% Input dirs
dset_name = 'ring1_zero_subset';

% Output dirs
output_name = '_indep_CG_TVphi1';
output_subdir = [dset_name,output_name];
num_ims = 10;
T = num_ims;

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);  

%Load first coupled fit
baseFileName = 'indep_fit_%i.mat';
load(fullfile(output_dir,sprintf(baseFileName,1)))
M = numel(dir(fullfile(output_dir,'*.mat')));

% Construct dictionary
switch P.basis
    case 'norm2'
        A0ft_stack = unshifted_basis_vector_ft_stack_norm2_zpad(P);
end
[N,K] = size(A0ft_stack);

%% Load data
B = zeros(N,T);

set = 1;
for j = 1:T
  b_data = load(fullfile(dataset,[P.prefix,'_',num2str(set),'.mat']));
    % Reduce image to vector if needed
    try
        b = sum(b_data.polar_image,1);
    catch
        b = b_data.polar_vector;
    end
    B(:,j) = b';
end

% Get error, sparsity, awmv
err_select = zeros(M,num_ims);
l0_select = zeros(M,num_ims);
l1_select = zeros(M,num_ims);
awmv_az = zeros(M,num_ims);
x_hats = cell(M,num_ims);
vdfs = zeros(K,T,M);
Fits = zeros(N,T,M);
vdf_time_all_fig = figure(56);
[ha3, pos3] = tight_subplot(1,M,[0.1 0.03],[.02 .08],[.02 .02]); 

for k = 1:M
    fprintf('%i of %i\n',k,M)
    % Load coupled fit 
    x_data = load(fullfile(output_dir,sprintf(baseFileName,k)));

    im_ind = 1;
    for j = 1:num_ims
        % Compute fit
        x = x_data.X_hat(:,:,j);
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x),P.params.zeroMask);
        Fits(:,j,k) = fit;
        
        % VDF/AWMV
        az_signal = squeeze(sum(x,1));
        var_sum = sum(az_signal(:));
        vdfs(:,j,k) = az_signal./var_sum;
        awmv_az(k,j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;

        % Objective terms
        err_select(k,j) = sum( (B(:,j)-fit(:)).^2 ) ;
        l0_select(k,j) = sum(x(:)>0);
        l1_select(k,j) = sum(x(:));

        im_ind = im_ind + 1;
    end
    tv_penalty(k) = sum(abs(DiffPhiX_1D(x_data.X_hat)),'all');

    % Plot VDFs
    axes(ha3(k))
    imagesc(squeeze(vdfs(:,:,k)))
    shading interp
    caxis([0 0.6])
    colormap(jet)

    ylabel('t')
    xlabel('\sigma')
end

%% Plot Fits
fits_fig = figure(222);
[ha2, pos2] = tight_subplot(2,5,[.005 .005],[.01 .01],[.01 .01]); 
ax_ind = 1;
for i = 1:num_ims
    axes(ha2(ax_ind))
    plot(B(:,i))
    hold on
    ax_ind = ax_ind + 1;
end
for k = 1:M
    ax_ind = 1;
    for i = 1:num_ims
        axes(ha2(ax_ind))
        plot(Fits(:,i,k))
        ax_ind = ax_ind + 1;
    end
end

%% Plot awmv
figure_dir = ['C:\Users\dpqb1\OneDrive\Desktop\',output_subdir,'_figures\'];
mkdir(figure_dir)

% Plot AWMV
awmv_fig = figure(1);
legend_str = {};

for k = 1:M
    hold on
    plot(awmv_az(k,:),'LineWidth',1.5)
    hold on

end
ylim([0 40])
ylabel('AWMV_\eta','FontSize',20)
xlabel('t','FontSize',20)
legend(legend_str,'location','best','FontSize',16)
saveas(awmv_fig,[figure_dir,'awmv_all_',dset_name,'.png'])
%% Save images
% 
% % Plot err
% total_err = sum(err_select,2);
% figure(2)
% hold on
% plot(0,sum(err_indep),'o')
% plot(lambda2_vals(sM:eM),total_err(sM:eM),'o-')
% ylabel('Error')
% xlabel('Coupling parameter')
% 
% % Plot l1 norm
% total_l1 = sum(l1_select,2);
% figure(3)
% hold on
% plot(0,sum(l1_indep),'o')
% plot(lambda2_vals(sM:eM),total_l1(sM:eM),'o-')
% ylabel('l1-norm')
% xlabel('Coupling parameter')
% 
% % Plot number nonzeros coefficients
% mean_l0 = mean(l0_select,2);
% figure(4)
% loglog(0,mean(l0_indep),'o')
% hold on
% plot(lambda2_vals(sM:eM),mean_l0(sM:eM),'o-')
% ylabel('l0-norm')
% xlabel('Coupling parameter')
% 
% 
% %     % Plot AWMV
% %     figure(5)
% %     legend_str = {};
% %     legend_str{1} = 'indep';
% %     hold on
% %     plot(awmv_az_init,'-','LineWidth',1.5)
% %     kk = 2;
% %     for k = [4,5,M]
% %         hold on
% %         plot(awmv_az(sort_i(k),:),'-','LineWidth',1.5)
% %         legend_str{kk} = sprintf('%0.04f',sort_lambda2(k));
% %         kk = kk + 1;
% %     end
% %     legend(legend_str,'Location','Best')
% 
% % % Plot L-curve 1
% % figure(6)
% % plot(obj2(sort_i),obj1(sort_i),'o-')
% % xlabel('l1-norm')
% % ylabel('Error')
% 
% % Plot L-curve 2
% 
% Lcurve_fig = figure(7);
% plot(total_err,tv_penalty,'o-')
% ylabel('TV')
% xlabel('Error')
% 
% total_err = sum(err_select(1:M,:),2);
% total_l1 = sum(tv_penalty(1:M,:),2);
% 
% % Find kink in L-cureve method #1
% slopes = (total_l1(2:end) - total_l1(1:end-1))./...
%      (total_err(2:end) - total_err(1:end-1));
% slope_select = find(abs(slopes)<1);
% select_ind = slope_select(1);
% select_ind = 26
% % Find kink in L-cureve method #2
% %     miny = min(tv_penalty);
% %     minx = min(total_err);
% %     coord = [total_err./minx,tv_penalty./miny];
% %     origin_dist = coord-1;
% %     [val,select_ind] = min(sum(origin_dist.^2,2));
% 
% selectx = total_err(select_ind);
% selecty = tv_penalty(select_ind,:);
% hold on
% plot(selectx,selecty,'s','Markersize',14)
% plot(sum(err_indep),tv_indep,'*','Markersize',14)
% saveas(Lcurve_fig,[figure_dir,'Lcurve_',dset_name,'_',dataset_num,'.png'])
% 
% % Plot paramter selected
% select_fig = figure(11);
% legend_str = {};
% legend_str{1} = 'truth';
% legend_str{2} = 'indep';
% hold on
% plot(awmv_truth,'LineWidth',1.5)
% plot(awmv_az_init,'LineWidth',1.5)
% kk = 3;
% for k = [17,select_ind,M]
%     hold on
%     plot(awmv_az(k,:),'LineWidth',1.5)
%     legend_str{kk} = sprintf('%0.01s',lambda2_vals(k));
%     kk = kk + 1;
% end
% ylim([0 40])
% ylabel('AWMV_\eta','FontSize',20)
% xlabel('t','FontSize',20)
% legend(legend_str,'location','best','FontSize',16)
% %     saveas(select_fig,[figure_dir,'awmv_select_',dset_name,'_',dataset_num,'.png'])
% 
% %% Plot vdf surface
% vdf_time_all_fig = figure(56);
% [ha3, pos3] = tight_subplot(3,ceil(M/3)+1,[0.1 0.03],[.02 .08],[.02 .02]); 
% im_ind = 1;
%     % Plot surface
% axes(ha3(im_ind))
% imagesc(squeeze(vdfs_indep'))
% shading interp
% caxis([0 0.6])
% colormap(jet)
% 
% title(['\lambda_2 = ','0'])
% ylabel('t')
% xlabel('\sigma')
% 
% im_ind = im_ind + 1;
% 
% for trial_k = sM:eM
%     fprintf('%i of %i\n',trial_k,M)
%     x_data = load(fullfile(output_dir,sprintf(baseFileName,trial_k)));
%     for image_num = 1:num_ims
%         x_hat = x_data.X_hat(:,:,image_num);
%         az_signal = squeeze(sum(x_hat,1));
%         var_sum = sum(az_signal(:));
%         vdf_time(trial_k,image_num,:) = az_signal/var_sum;
%     end
% 
%     % Plot surface
%     axes(ha3(im_ind))
%     imagesc(squeeze(vdf_time(trial_k,:,:)))
%     shading interp
%     caxis([0 0.6])
%     colormap(jet)
% 
%     title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(trial_k))])
%     ylabel('t')
%     xlabel('\sigma')
% 
%     im_ind = im_ind + 1;
% end
% saveas(vdf_time_all_fig,[figure_dir,'vdf_time_all_',dset_name,'_',dataset_num,'.png'])
% 
% vdf_time_fig = figure(566);
% 
% % Plot surface
% imagesc(squeeze(vdf_time(select_ind,:,:)))
% shading interp
% caxis([0 0.6])
% colormap(jet)
% 
% title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(select_ind))])
% ylabel('t')
% xlabel('\sigma')
% %     saveas(vdf_time_fig,[figure_dir,'vdf_time_select_',dset_name,'_',dataset_num,'.png'])
% 
% vdf_indep_fig = figure(567);
% % Plot surface
%     imagesc(squeeze(vdfs_indep'))
% shading interp
% caxis([0 0.6])
% colormap(jet)
% 
% title(['\lambda_2 = ','0'])
% ylabel('t')
% xlabel('\sigma')
