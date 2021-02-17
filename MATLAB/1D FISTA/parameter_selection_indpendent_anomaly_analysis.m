clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;
% Parent directory
top_dir = 'D:\CHESS_data';
%     top_dir = '/cluster/shared/dbanco02';

% Input dirs
for data_num = 2
    close all
dset_name = ['simulated_two_spot_1D_anomaly_',num2str(data_num)];
noise_added = 0:0.03:0.30;
noise_thresh = [0.01,0.03:0.03:0.24,0.3,0.5];

% Output dirs
output_name = '_indep_ISM1';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);

baseFileName = 'indep_fit_%i_%i.mat';

% Load most parameters by loading single output
load(fullfile(output_dir,sprintf(baseFileName,1,1)))

% Real lambda values
lambda_values = P.lambda_values;
% lambda_values = [logspace(-7,-5.1,10) logspace(-5,1,30)];
% P.lambda_values = lambda_values;

N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;
M = numel(lambda_values);

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

nosie_est = norm(randn(179,1)*0.03).^2;
for i = 1:M
    fprintf('%i of %i \n',i,M)
    for j = 1:T
        b_data = load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']));
        polar_vector = b_data.polar_vector(1:179);
        e_data = load(fullfile(output_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        fit = Ax_ft_1D(A0ft_stack,e_data.x_hat);
        b = P.dataScale*polar_vector;
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

%% L curve parameter selection for l1-norm term
select_indices = zeros(T,1);
for t = 1:T
    err_t = err_select(:,t);
    l1_t = l1_select(:,t);
    err_t = err_t;%/max(err_t);
    l1_t = l1_t;%/max(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist + (err_t == 0) )  );
end

for t = 1:T
    load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    b = P.dataScale*polar_vector(1:179);
    rel_err_t = err_select(:,t)/sum(b(:).^2);
    while rel_err_t(select_indices(t)) > noise_thresh(data_num)
        if select_indices(t) > 1
            select_indices(t) = select_indices(t) - 1;
        else
            select_indices(t) = find(rel_err_t == min(rel_err_t));
            break
        end
    end
end


% Selected VDF;
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

% Plot AWMV

% Load True AWMV
true_name = 'simulated_two_spot_1D_anomaly_11';
true_data = load(fullfile(top_dir,true_name,'synth_data.mat'));
true_awmv = zeros(T,1);
for t = 1:T
    true_stds = true_data.synth_sample{t}.std_theta;
    true_amps = true_data.synth_sample{t}.amplitudes;
    true_awmv(t) = true_stds'*true_amps/sum(true_amps);
end


awmv_fig = figure(9);
hold on
plot(true_awmv)
plot(select_vdf*sqrt(P.var_theta)')
shading interp
caxis([0 0.6])
colormap(jet)
title('VDF selected parameters')
legend('Truth','Indep','Location','Best')
saveas(awmv_fig,['C:\Users\dpqb1\Desktop\indep_param_select\awmv_fig_',num2str(data_num),'.png'])

% Plot fits for L-curve selected parameter
fits_fig = figure(10);
[ha2, ~] = tight_subplot(4,5,[.005 .005],[.01 .01],[.01 .01]); 
awmv_az_vdfs = zeros(T,1);
im_ind = 1;
for t = 1:T
    load(fullfile(output_dir,sprintf(baseFileName,select_indices(t),t)))
    load(fullfile(dataset,[P.prefix,'_',num2str(t),'.mat']))
    
    % remove small coefs
    final_thresh = 1e-3*sum(x_hat(:));
    x_hat(x_hat<final_thresh) = 0;
   
    polar_vector = polar_vector(1:179);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_az_vdfs(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    b = P.dataScale*zeroPad(polar_vector,P.params.zeroPad);
    

    % Plot
    axes(ha2(im_ind))
    hold on
    plot(b)
    plot(fit)
    rel_err = sum((fit(:)-b(:)).^2)./sum(b(:).^2);
    legend(sprintf('%i \n %0.2f',sum(x_hat(:)>0),rel_err),'location','northeast')
    im_ind = im_ind + 1;
end
saveas(fits_fig,['C:\Users\dpqb1\Desktop\indep_param_select\fits_fig_',num2str(data_num),'.png'])
end
%% Show coefficients of selected parameters
% fits_fig = figure(11);
% [ha2, ~] = tight_subplot(2,T/2,[.005 .005],[.01 .01],[.01 .01]); 
% awmv_az_vdfs = zeros(T,1);
% im_ind = 1;
% for t = 1:T
%     load(fullfile(output_dir,sprintf(baseFileName,select_indices(t),t)))
%     axes(ha2(t))
%     imagesc(x_hat)
% end

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