%% Parameter selection
clear all
close all

% Parameter selection
disp('Setup parms')
P.set = 1;
% Parent directory
top_dir = 'D:\CHESS_data';

% Input dirs
for dd_num = 5
% dd_num = 2;
    close all
dset_name = ['simulated_two_spot_1D_anomaly_',num2str(dd_num)];
% figure_dir = ['C:\Users\dpqb1\Desktop\anomaly_figures_2norm3\'];
figure_dir = ['C:\Users\dpqb1\Desktop\anomaly_figures_3\'];
mkdir(figure_dir)
% Indep dirs
indep_name = '_indep_ISM1';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
% output_name = '_coupled_CG_TVphi_2norm3';
output_name = '_coupled_CG_TVphi3'
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);

% File Parameters
baseFileName = 'coupled_fit_%i.mat';

% Universal Parameters
% Ring sampling parameters
load(fullfile(output_dir,sprintf(baseFileName,1)))
lambda2_vals = P.lambda2_values(1:end);
% lambda2_vals = [logspace(-6,-4.1,15) logspace(-4,1,30)]
M = numel(lambda2_vals);

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

%%
T = P.num_ims;
[N,K] = size(A0ft_stack);
B = zeros(N,T);

for k = 1:M
    fprintf('%i of %i\n',k,M)
    x_data = load( fullfile(output_dir,sprintf(baseFileName,k)) );
    
    for j = 1:P.num_ims
        % Fit objective
        x = x_data.X_hat(:,:,j);
        fit = Ax_ft_1D(A0ft_stack,x);
        
        if k == 1
            load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']) )
            b = P.dataScale*polar_vector(1:179);
            B(:,j) = b';
        else
            b = B(:,j);
        end
        err_select(k,j) = sum( (b(:)-fit(:)).^2 )./sum( b(:).^2 ) ;

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
    x = e_data.x_hat;
    X_indep(:,:,j) = x;
    fit = Ax_ft_1D(A0ft_stack,x);
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

% Load True AWMV
true_name = 'simulated_two_spot_1D_anomaly_11';
true_data = load(fullfile(top_dir,true_name,'synth_data.mat'));
true_awmv = zeros(T,1);
for t = 1:T
    true_stds = true_data.synth_sample{t}.std_theta;
    true_amps = true_data.synth_sample{t}.amplitudes;
    true_awmv(t) = true_stds'*true_amps/sum(true_amps);
end

%% Plot awmv

% Plot AWMV
all_awmv_fig = figure(1);
legend_str = {};
% legend_str{1} = 'truth';
legend_str{1} = '0';
kk = 2;
hold on
plot(true_awmv,'LineWidth',1.5)
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

total_err = sum(err_select,2);
total_l1 = sum(l1_select,2);
mean_l0 = mean(l0_select,2);

saveas(all_awmv_fig,[figure_dir,'all_awmv',dset_name,'_',num2str(dd_num),'.png'])

%% Save images

Lcurve_fig = figure(7);
loglog(total_err,tv_penalty,'o-')
ylabel('TV')
xlabel('Error')

total_err = sum(err_select,2);

% Find kink in L-cureve method #3
err_scale = total_err/max(total_err(:));
tv_scale = tv_penalty/max(tv_penalty(:));
sq_origin_dist = abs(tv_scale).^2 + abs(err_scale).^2;
select_ind = find(sq_origin_dist == min(sq_origin_dist),1);

% Just use closest to truth
% awmv_error = sum((awmv_az-repmat(true_awmv,[1,30])').^2,2);
% min_err_ind = find(awmv_error == min(awmv_error));
% select_ind = min_err_ind(1);

selectx = total_err(select_ind);
selecty = tv_penalty(select_ind,:);
hold on
plot(selectx,selecty,'s','Markersize',14)
plot(sum(err_indep),tv_indep,'*','Markersize',14)

saveas(Lcurve_fig,[figure_dir,'L_curve',dset_name,'_',num2str(dd_num),'.png'])


figure(8)
plot(lambda2_vals,tv_penalty,'o-')
ylabel('tv')
xlabel('Coupling parameter')

%% Plot paramter selected
select_fig = figure(11);
legend_str = {};
legend_str{1} = 'truth';
legend_str{2} = 'indep';
hold on
plot(true_awmv,'LineWidth',1.5)
plot(awmv_az_init,'LineWidth',1.5)
kk = 3;
%  select_ind = 1;
for k = [select_ind]
    hold on
    
    plot(awmv_az(k,:),'LineWidth',1.5)
    legend_str{kk} = sprintf('%0.01s',lambda2_vals(k));
    kk = kk + 1;
end

ylabel('AWMV_\eta','FontSize',20)
xlabel('t','FontSize',20)
legend(legend_str,'location','best','FontSize',16)
saveas(select_fig,[figure_dir,'awmv_select_',dset_name,'_',num2str(dd_num),'.png'])
save([dset_name,'_couple_fit_data.mat'],'awmv_az_init','awmv_az','select_ind')


%% Plot fits
fits_fig = figure(222);
[ha2, pos2] = tight_subplot(5,4,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
x_data = load(fullfile(output_dir,sprintf(baseFileName,select_ind)));
for image_num = 1:T
    x_hat = x_data.X_hat(:,:,image_num);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = B(:,image_num);
    
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',image_num),'location','northeast')
    im_ind = im_ind + 1;
end
saveas(fits_fig,[figure_dir,'fits_',dset_name,'_',num2str(dd_num),'.png'])

%% Plot vdfs and indep vdfs
fits_fig = figure(333);
[ha2, pos2] = tight_subplot(5,4,[.005 .005],[.01 .01],[.01 .01]); 
unnorm_vdfs = zeros(K,T);
unnorm_vdfs_indep = zeros(K,T);
im_ind = 1;
x_data = load(fullfile(output_dir,sprintf(baseFileName,select_ind)));
for image_num = 1:T
    x_hat = x_data.X_hat(:,:,image_num);
    az_signal = squeeze(sum(x_hat,1));
    unnorm_vdfs(:,image_num) = az_signal(:);
    x_hat = X_indep(:,:,image_num);
    az_signal = squeeze(sum(x_hat,1));
    unnorm_vdfs_indep(:,image_num) = az_signal(:);
end
for k = 1:K
    axes(ha2(k))
    plot(unnorm_vdfs(k,:),'-o')
    hold on
    plot(unnorm_vdfs_indep(k,:),'-o')
end
% saveas(fits_fig,[figure_dir,'fits_',dset_name,'_',num2str(dd_num),'.png'])

figure(45)
imagesc(unnorm_vdfs)
xlabel('time')
ylabel('\sigma_\eta')

figure(46)
imagesc(unnorm_vdfs_indep)
xlabel('time')
ylabel('\sigma_\eta')

indep_diff1 = 0;
indep_diff2 = 0;
diff1 = 0;
diff2 = 0;
for t = 1:T-1
    diff1 = diff1 + sum(abs(unnorm_vdfs(:,t+1)-unnorm_vdfs(:,t)));
    indep_diff1 = indep_diff1 + sum(abs(unnorm_vdfs_indep(:,t+1)-unnorm_vdfs_indep(:,t)));
    diff2 = diff2 + sum((unnorm_vdfs(:,t+1)-unnorm_vdfs(:,t)).^2);
    indep_diff2 = indep_diff2 + sum((unnorm_vdfs_indep(:,t+1)-unnorm_vdfs_indep(:,t)).^2);
end
diff1
diff2
indep_diff1
indep_diff2

end

