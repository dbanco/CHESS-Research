%% Parameter selection
clear all
close all

% Parameter selection
disp('Setup parms')
P.set = 1;
% Parent directory
top_dir = 'E:\PureTiRD_full';
% top_dir = 'E:\MMPAD_data';
% Input dirs
dset_name = 'ring2_zero';

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
% lambda2_vals = [logspace(-6,-4.1,15) logspace(-4,1,30)]
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
        if size(polar_image,1) == 261
            aaa = polar_image';
            clear polar_image
            polar_image = aaa;
            save(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']),'polar_image')
        end
        b = P.dataScale*sum(polar_image,1);
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

%% Compute tv objective
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
    b = P.dataScale*sum(polar_image,1);
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
legend_str{1} = '0';
kk = 2;
hold on
plot(awmv_az_init,'LineWidth',1.5)
for k = 19
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


% Find kink in L-cureve method #1
% slopes = (total_l1(2:end) - total_l1(1:end-1))./...
%      (total_err(2:end) - total_err(1:end-1));
% slope_select = find(abs(slopes)<1);
% select_ind = slope_select(1);

% Find kink in L-cureve method #2
%     miny = min(tv_penalty); minx = min(total_err); coord =
%     [total_err./minx,tv_penalty./miny]; origin_dist = coord-1;
%     [val,select_ind] = min(sum(origin_dist.^2,2));

% Find kink in L-cureve method #3
err_scale = total_err/max(total_err(:));
tv_scale = tv_penalty/max(tv_penalty(:));
sq_origin_dist = abs(tv_scale).^2 + abs(err_scale).^2;
select_ind = find(sq_origin_dist == min(sq_origin_dist),1);
select_ind = 22;
selectx = total_err(select_ind);
selecty = tv_penalty(select_ind,:);
hold on
plot(selectx,selecty,'s','Markersize',14)
plot(sum(err_indep),tv_indep,'*','Markersize',14)


figure(8)
plot(lambda2_vals,tv_penalty,'o-')
ylabel('tv')
xlabel('Coupling parameter')


% figure(77)
% plot(err_scale,tv_scale,'o-')
% ylabel('TV')
% xlabel('Error')

figure(9)
plot(err_select(select_ind,:))
hold on
plot(err_indep)
legend('coupled','indep')
ylabel('error')
xlabel('time')

figure(12) % analyze correspondence of l1-parameters to awmv
cutoff = mean(awmv_az_init(200:end));
param_cut = (P.params.lambda1/max(P.params.lambda1(:))*max(awmv_az_init)) > cutoff;
indep_cut = awmv_az_init > cutoff;
correspond = param_cut(:) == indep_cut(:);
% plot(P.params.lambda1/max(P.params.lambda1(:))*max(awmv_az_init))
% plot(200:546,awmv_az_init(200:end).*correspond(200:end),'o')
hold on
plot(awmv_az_init)
plot(awmv_az(select_ind,:))
legend('indep','coupled','Location','Best')

figure(122) % analyze correspondence of l1-parameters to awmv
cutoff = mean(awmv_az_init(200:end));
param_cut = (P.params.lambda1/max(P.params.lambda1(:))*max(awmv_az_init)) > cutoff;
indep_cut = awmv_az_init > cutoff;
correspond = param_cut(:) == indep_cut(:);
% plot(P.params.lambda1/max(P.params.lambda1(:))*max(awmv_az_init))
% plot(200:546,awmv_az_init(200:end).*correspond(200:end),'o')
hold on
plot(awmv_az_init)
plot(awmv_az(select_ind,:))
title('All scaled')
legend('parameter value','indep','coupled','Location','Best')

figure(13) % analyze correspondence of l1-parameters to awmv
indep_end = awmv_az_init(200:end);
coupled_end = awmv_az(select_ind,200:end);
params_end = P.params.lambda1(200:end)*20 + mean(indep_end);

cutoff = mean(awmv_az_init(200:end));
param_cut = (P.params.lambda1/max(P.params.lambda1(:))*max(awmv_az_init)) > cutoff;
indep_cut = awmv_az_init > cutoff;
correspond = param_cut(:) == indep_cut(:);
plot(params_end)
% plot(200:546,awmv_az_init(200:end).*correspond(200:end),'o')
hold on
plot(indep_end)
plot(coupled_end)
title('All scaled')
legend('parameter value','indep','coupled','Location','Best')

%% Plot paramter selected
select_fig = figure(11);
legend_str = {};
legend_str{1} = 'indep';
hold on
plot(awmv_az_init,'LineWidth',1.5)
kk = 2;

for k = [select_ind]
    hold on
    plot(awmv_az(k,:),'LineWidth',1.5)
    legend_str{kk} = sprintf('%0.01s',lambda2_vals(k));
    kk = kk + 1;
end

ylabel('AWMV_\eta','FontSize',20)
xlabel('t','FontSize',20)
legend(legend_str,'location','best','FontSize',16)
%     saveas(select_fig,[figure_dir,'awmv_select_',dset_name,'_',dataset_num,'.png'])
save([dset_name,'_couple_fit_data.mat'],'awmv_az_init','awmv_az','select_ind')
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
% saveas(vdf_time_all_fig,[figure_dir,'vdf_time_all_',dset_name,'_',dataset_num,'.png'])

%% 
vdf_time_fig = figure(566);

subplot(2,1,1)
% Plot surface
imagesc(squeeze(vdf_time(select_ind,:,:))')
shading interp
caxis([0 0.6])
colormap(jet)
colorbar()
title(['\lambda_2 = ',sprintf('%1.1d',lambda2_vals(select_ind))])
xlabel('t')
ylabel('\sigma')
%     saveas(vdf_time_fig,[figure_dir,'vdf_time_select_',dset_name,'_',dataset_num,'.png'])

subplot(2,1,2)
% Plot surface
    imagesc(squeeze(vdfs_indep))
shading interp
caxis([0 0.6])
colormap(jet)
colorbar()
title(['\lambda_2 = ','0'])
xlabel('t')
ylabel('\sigma')

%% Plot indep fits
fits_fig = figure(2222);
[ha2, pos2] = tight_subplot(10,7,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
% x_data = load(fullfile(output_dir,sprintf(baseFileName,select_ind)));
for image_num = 1:T
    x_hat = X_indep(:,:,image_num);
    fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_hat),P.params.zeroMask);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(image_num) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    load(fullfile(dataset,[P.prefix,'_',num2str(image_num),'.mat']) )
    polar_vector = sum(polar_image,1);
    b = P.dataScale*zeroPad(polar_vector,P.params.zeroPad);
    
    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    legend(sprintf('%i',image_num),'location','northeast')
    im_ind = im_ind + 1;
end


%% Plot fits
fits_fig = figure(222);
[ha2, pos2] = tight_subplot(10,10,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
x_data = load(fullfile(output_dir,sprintf(baseFileName,select_ind)));
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
    legend(sprintf('%i',image_num),'location','northeast')
    im_ind = im_ind + 1;
end
%     saveas(fits_fig,[figure_dir,'fits_',dset_name,'_',dataset_num,'.png'])

%% Try connected components analysis
image_num = 1;
thresh = 1e-3;
x_hat = x_data.X_hat(:,:,image_num);
shift_x = shift2D(x_hat,130,0);
cc = bwconncomp(shift_x>thresh);
cc_reg = regionprops(cc,'Centroid');

colors = jet(cc.NumObjects);

load(fullfile(dataset,[P.prefix,'_',num2str(image_num),'.mat']) )
polar_vector = sum(polar_image,1);
b = P.dataScale*zeroPad(polar_vector,P.params.zeroPad);

figure(6)
hold on
plot(b)
fit = Ax_ft_1D(A0ft_stack,x_hat);
plot(fit)

for i = 1:cc.NumObjects
    aa = cc_reg(i).Centroid(2);
    bb = cc_reg(i).Centroid(2);
    patch = b(floor(aa):ceil(bb));
    patch_max = max(patch(:));
    plot([aa,bb],[1,1]*patch_max,'o-','Color',colors(i,:))

end

% Compute AWMV for each connected component
awmv_peaks = zeros(cc.NumObjects,1);
for i = 1:cc.NumObjects
    pixelList = cc.PixelIdxList{i};
    for j = 1:numel(pixelList)
        [row,col] = ind2sub(size(shift_x),pixelList(j));
        awmv_peaks(i) = awmv_peaks(i) + shift_x(row,col)*sqrt(P.var_theta(col));
    end
    aa = cc_reg(i).Centroid(2);
    bb = cc_reg(i).Centroid(2);
    patch = b(floor(aa):ceil(bb));
    patch_max = max(patch(:));
    plot([aa,bb],[1,1]*patch_max,'o-','Color',colors(i,:))

end

figure(8)
subplot(2,1,1)
imagesc(shift_x>thresh)
subplot(2,1,2)
bar(awmv_peaks)

%% Evaluate awmv of individual peaks in time
awmv_peaks = zeros(15,200);
for image_num = 1:200
    thresh = 1e-3;
    x_hat = x_data.X_hat(:,:,image_num);
    shift_x = shift2D(x_hat,130,0);
    cc = bwconncomp(shift_x>thresh);
    cc_reg = regionprops(cc,'Centroid');
    
    % Compute AWMV for each connected component
    for i = 1:cc.NumObjects
        pixelList = cc.PixelIdxList{i};
        for j = 1:numel(pixelList)
            [row,col] = ind2sub(size(shift_x),pixelList(j));
            awmv_peaks(i,image_num) = awmv_peaks(i) + shift_x(row,col)*sqrt(P.var_theta(col));
        end
    end
end

figure(6)
imagesc(awmv_peaks)
colorbar()

%% Track particular connected component (doesnt track same location?)
figure(7)
cc_num = 7;
b_single_cc = zeros(numel(b),200);
thresh = 1e-3;
centroid1 = 0;
for image_num = 1:200
    
    x_hat = x_data.X_hat(:,:,image_num);
    shift_x = shift2D(x_hat,130,0);
    cc = bwconncomp(shift_x>thresh);
    cc_reg = regionprops(cc,'Centroid');
    
    % Initialize centroid
    if image_num == 1
        centroid1 = cc_reg(cc_num).Centroid(2);
    end
    centroids = reshape([cc_reg.Centroid],[2,cc.NumObjects]);
    centroid_locs = centroids(2,:);
    [~,cc_idx] = min(abs(centroid1-centroid_locs));
    
    % Update centroid
    if image_num ~= 1
        centroid1 = centroid_locs(cc_idx);
    end
    
    % Plot selected connected component
    x_single_cc = zeros(size(x_hat));
    pixelList = cc.PixelIdxList{cc_idx};    
    x_single_cc(pixelList) = x_hat(pixelList);
    b_single_cc(:,image_num) = Ax_ft_1D(A0ft_stack,x_single_cc);
    plot(b_single_cc(:,image_num))
    pause(0.1)
end

%{

%% Analyze some stuff
fits_fig = figure(888);
[ha8, pos8] = tight_subplot(10,10,[.005 .005],[.01 .01],[.01 .01]); 
im_ind = 1;
subT = 100;
x_data = load( fullfile(output_dir,sprintf(baseFileName,select_ind)) );
err_b = zeros(subT,1);
l1_b = zeros(subT,1);
l0_b = zeros(subT,1);
bmean = zeros(subT,1);
awmv_az_b = zeros(subT,1);
vdfs_b = zeros(K,subT);
cutoff = K;
for j = 201:300
    % Fit objective
    x = x_data.X_hat(:,:,j);
    x1 = x;
    x2 = x;
    x1(:,cutoff) = 0;
    x2(:,1:(cutoff-1)) = 0;
%     fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x),P.params.zeroMask);
%     fit1 = forceMaskToZero(Ax_ft_1D(A0ft_stack,x1),P.params.zeroMask);
%     fit2 = forceMaskToZero(Ax_ft_1D(A0ft_stack,x2),P.params.zeroMask);
    
    load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']) )
    b = P.dataScale*sum(polar_image,1);
    err_b(im_ind) = sum( (b(:)-fit(:)).^2 )./sum( b(:).^2 ) ;
    bnorms(im_ind) = norm(b);
    bmean(im_ind) = mean(b);
    % Vdf
    az_signal = squeeze(sum(x,1));
    var_sum = sum(az_signal(:));
    vdfs_b(:,im_ind) = az_signal./var_sum;
    awmv_az_b(im_ind) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;

    % Sparsity objective
    l0_b(im_ind) = sum(x(:)>0);
    l1_b(im_ind) = sum(x(:));
    
    figure(8989)
%     axes(ha8(im_ind))
    hold on
%     plot(b)
%     plot(fit1)
    cmap = hot(100);
    plot(b,'Color',cmap(im_ind,:))
    legend(sprintf('%i',j),'location','northeast')
    ylim([0 0.25])
    im_ind = im_ind + 1;
    
end

num_plots = 6;
figure(88)
subplot(num_plots,1,1)
plot(err_select(select_ind,:))
title('error')
subplot(num_plots,1,2)
plot(l1_b)
title('l1-norm')
subplot(num_plots,1,3)
plot(bnorms)
title('||b||_2')
subplot(num_plots,1,4)
plot(P.params.lambda1)
title('\lambda_1')
subplot(num_plots,1,5)
plot(bmean)
title('bmean')
subplot(num_plots,1,6)
plot(awmv_az_b)
title('awmv')
figure(89)
imagesc(vdfs_b,[0 0.05])
%}
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

