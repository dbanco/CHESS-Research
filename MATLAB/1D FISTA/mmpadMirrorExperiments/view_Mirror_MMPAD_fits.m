clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;
% Parent directory
top_dir = 'D:\MMPAD_data_nr1';
% top_dir = 'E:\PureTiRD_full';
%     top_dir = '/cluster/shared/dbanco02';

for ring_num = 3
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
%{
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
fits_indep = zeros(N,T);
for j = 1:T
    i_data = load(fullfile(top_dir,indepDir,sprintf(indepName,e_data.P.params.lambda1_indices(j),j)),'err','x_hat');
    load(fullfile(dataset,[P.prefix,'_',num2str(j),'.mat']) )
    x = i_data.x_hat;
    X_indep(:,:,j) = x;
    fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,x),P.params.zeroMask);
    fits_indep(:,j) = fit;
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
%}
%% Plot fits for L-curve selected parameter
load(fullfile(output_dir,sprintf(baseFileName,25)))
fits_fig = figure(99);
plot_rows = 12;
plot_cols = 5; 
start_t = 1;
end_t = start_t + plot_rows*plot_cols - 1;
% [ha2, ~] = tight_subplot(plot_rows,plot_cols,[.005 .005],[.01 .01],[.01 .01]); 
im_ind = 1;
T_end = 201;
fits_coupled = zeros(N,T_end);
for t = 1:T_end
    x_hat = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    fits_coupled(:,t) = fit;
%     az_signal = squeeze(sum(x_hat,1));
%     var_sum = sum(az_signal(:));
%     b = B(:,t);

%     % Plot
%     axes(ha2(im_ind))
%     hold on
%     plot(b((pad1+1):(pad1+nn)))
%     plot(fit((pad1+1):(pad1+nn)))
%     rel_err = sum((fit(:)-b(:)).^2)/norm(b(:))^2;
%     legend(sprintf('%0.2f',(t-1)/4),'location','best')
%     im_ind = im_ind + 1;
end

% MMPAD metadata
mmpad_dims = [396,265];
rings = {'020','112','021','004'};
d_space = [1.27773,1.24845,1.23272,1.17138];
two_theta = [6.89132,7.05316,7.14328,7.5179];
% eta axis
pixel_side = 150e-6;
gap_width = 0.75e-3;
detec_dist = 4.6;
detect_angle = 14.6;
% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;
% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;
% time axis
total_time = 136.25;
fps = 4;
x_time = 0:1/fps:total_time;
eta_range = linspace(-261/2,261/2,261)*pixel_angle(1) + detect_angle;
N_eta = numel(eta_range);
clim1 = 0;
clim2 = 5;

% Data water fall for first 25s
figure(1)
b_z = (B(pad1+1:N-pad2,1:T_end))';
h = surf((eta_range),(x_time(1:T_end)),b_z);

ylabel('time(s)','FontSize',16)
xlabel('\eta (\circ)','FontSize',18)
zlabel('log(Intensity)','FontSize',16)
set(gca, 'ZScale', 'log');
set(gca, 'ColorScale', 'log');
set(h,'linestyle','none');
zlim([clim1,clim2])
set(gca,'CLim', [clim1 clim2])
hold on
z_max = max(b_z(:));

line(eta_range,ones(N_eta,1)*10,log(z_max)*ones(N_eta,1),'Linewidth',1,'Color','r')
% set(h,'linestyle','none');
% % Coupled waterfall for first 25s
% figure(2)
% waterfall((eta_range),(x_time(1:T_end)),(fits_coupled(pad1+1:N-pad2,1:T_end)'))
% ylabel('time(s)','FontSize',16)
% xlabel('\eta (\circ)','FontSize',18)
% zlabel('log(Intensity)','FontSize',16)
% set(gca, 'ZScale', 'log');
% set(gca, 'ColorScale', 'log');
% zlim([clim1,clim2])
% set(gca,'CLim', [clim1 clim2])
% 
% % Indep waterfall for first 25s
% figure(3)
% waterfall((eta_range),(x_time(1:T_end)),(fits_indep(pad1+1:N-pad2,1:T_end)'))
% ylabel('time(s)','FontSize',16)
% xlabel('\eta (\circ)','FontSize',18)
% zlabel('log(Intensity)','FontSize',16)
% set(gca, 'ZScale', 'log');
% set(gca, 'ColorScale', 'log');
% zlim([clim1,clim2])
% set(gca,'CLim', [clim1 clim2])

end