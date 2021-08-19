clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;
% Parent directory
top_dir = 'D:\MMPAD_data_nr1';
% top_dir = 'E:\PureTiRD_full';
%     top_dir = '/cluster/shared/dbanco02';

for ring_num = 4
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

%% Plot fits for L-curve selected parameter
load(fullfile(output_dir,sprintf(baseFileName,25)))
fits_fig = figure(99);
plot_rows = 12;
plot_cols = 5; 
start_t = 1;
end_t = start_t + plot_rows*plot_cols - 1;
% [ha2, ~] = tight_subplot(plot_rows,plot_cols,[.005 .005],[.01 .01],[.01 .01]); 
im_ind = 1;

line_t1 = 25;
line_t2 = 50;
line_t3 = 36;
T_end = 201;
t_range = (line_t1*4+1):(line_t2*4+1);


fits_coupled = zeros(N,T_end);
for t = 1:T_end
    x_hat = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    fits_coupled(:,t) = fit;
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
clim1 = 2e-1;
clim2 = 1e3;

% Data water fall for first 25s
figure(1)
b_z = (B(pad1+1:N-pad2,1:T_end))';
h = surf((eta_range),(x_time(1:T_end)),b_z);
shading flat
colorbar()
ylabel('time(s)','FontSize',16)
xlabel('\eta (\circ)','FontSize',18)
zlabel('log(Intensity)','FontSize',16)
set(gca, 'ZScale', 'log');
set(gca, 'ColorScale', 'log');
set(h,'linestyle','none');

% zlim([clim1,clim2])
% set(gca,'CLim', [clim1 clim2])
hold on
z_max = max(b_z(:));

line(eta_range,ones(N_eta,1)*line_t1,log(z_max)*ones(N_eta,1),'Linewidth',1,'Color','r')
line(eta_range,ones(N_eta,1)*line_t2,log(z_max)*ones(N_eta,1),'Linewidth',1,'Color','r')
line(eta_range,ones(N_eta,1)*line_t3,z_max*ones(N_eta,1),'Linewidth',1,'Color','r')

figure(2)
clim1 = 0;
clim2 = 2;
b_z = (B(pad1+1:N-pad2,t_range))';
hh = surf((eta_range),(x_time(t_range)),b_z);
shading flat
colorbar()
ylabel('time(s)','FontSize',16)
xlabel('\eta (\circ)','FontSize',18)
zlabel('log(Intensity)','FontSize',16)
% set(gca, 'ZScale', 'log');
% set(gca, 'ColorScale', 'log');
set(hh,'linestyle','none');

% zlim([clim1,clim2])
% set(gca,'CLim', [clim1 clim2])
hold on
z_max = max(b_z(:));

line(eta_range,ones(N_eta,1)*line_t1,z_max*ones(N_eta,1),'Linewidth',1,'Color','r')
line(eta_range,ones(N_eta,1)*line_t2,z_max*ones(N_eta,1),'Linewidth',1,'Color','r')
line(eta_range,ones(N_eta,1)*line_t3,z_max*ones(N_eta,1),'Linewidth',1,'Color','r')

%% Plot all coupled awmv
datadir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
fileBase = ['ring%i_zero_mirror_coupled_awmv.mat'];
outName =  ['allAWMV_nr1.png'];
fig_out = figure(5);
R = 4;
select_inds = [29 23 25 25];
for ring = 1:R

    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(x_time,awmv_az(select_inds(ring),:)*pixel_angle(ring),'LineWidth',1.5)
    ylabel('AWMV in \eta (degrees)')
    xlabel('time (s)')

end
c_leg = legend(rings,'Location','Best');
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')
grid on

line(ones(N_eta,1)*line_t1,linspace(0,1.3,N_eta),'Linewidth',1,'Color','r')
line(ones(N_eta,1)*line_t2,linspace(0,1.3,N_eta),'Linewidth',1,'Color','r')
line(ones(N_eta,1)*line_t3,linspace(0,1.3,N_eta),'Linewidth',1,'Color','r')

saveas(fig_out,fullfile(datadir,outName))
end