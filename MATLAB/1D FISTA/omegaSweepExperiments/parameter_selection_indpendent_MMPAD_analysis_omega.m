clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;

% Parent directory
% top_dir = 'E:\MMPAD_omega';
top_dir = '/cluster/shared/dbanco02/data/MMPAD_omega';
om_dir = {'omega2','omega3','omega4','omega5'};
r_dir = {'ring1','ring2','ring3','ring4'};

for o =1:4
    for r = 1:4
        close all

% Input dirs
dset_name = r_dir{r};
om_name = om_dir{o};

% Output dirs
output_name = '_indep_ISM_Mirror4';
output_subdir = [dset_name,om_name,output_name];

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

[N,K] = size(x_hat);
% T = P.num_ims;
T = 546;
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
rel_err_select = zeros(M,T);
l0_select = zeros(M,T);
l1_select = zeros(M,T);
x_indep = cell(T,1);
awmv = zeros(T,M);

B = zeros(N,T);
% Load data
for j = 1:T
    b_data = load(fullfile(top_dir,om_name,dset_name,[P.prefix,'_',num2str(j),'.mat']));
    b = P.dataScale*sum(b_data.polar_image,2);
    
    % Mirror data
    nn = numel(b);
    pad1 = floor(nn/2);
    pad2 = ceil(nn/2);
    N = nn + pad1 + pad2;
    b_mirror = zeros(N,1);
    b_mirror((pad1+1):(pad1+nn)) = b;
    b_mirror((1+N-pad2):N) = flip(b((nn-pad2+1):nn));
    b_mirror(1:pad1) = flip(b(1:pad1));
    B(:,j) = b_mirror;
end

tv_time = zeros(M,T-1);
im_ind = 1;
for i = 1:M
    fprintf('%i of %i \n',i,M)
    for j = 1:T
        b = B(:,j);
        e_data = load(fullfile(output_dir,sprintf(baseFileName,i,j)),'err','x_hat');
        e_data.x_hat(e_data.x_hat<0) = 0;
        
        fit = Ax_ft_1D(A0ft_stack,e_data.x_hat);   
        N = size(e_data.x_hat,1);
        err_select(i,j) = sum(( fit(:) - b(:) ).^2);
        rel_err_select(i,j) = sqrt(err_select(i,j))/norm(b(:));
        l0_select(i,j) = sum(e_data.x_hat(:) > 1e-4*sum(e_data.x_hat(:)));
        l1_select(i,j) = sum(abs(e_data.x_hat(:)));
        az_signal = squeeze(sum(e_data.x_hat,1));
        var_sum = sum(az_signal(:));
        vdf_time(i,j,:) = az_signal/var_sum;
        awmv(j,i) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        
    end

    im_ind = im_ind + 1;
end
err_select(err_select > 10^10) = 0;
l0_select(l0_select > 10^10) = 0;
l1_select(l1_select > 10^10) = 0;

%% L curve parameter selection
select_indices = zeros(T,1);
selected_error = zeros(T,1);
for t = 1:T
    err_t = err_select(1:M,t);
    l1_t = l1_select(1:M,t);
    err_t = err_t/max(err_t);
    l1_t = l1_t/max(l1_t);
    sq_origin_dist = abs(l1_t) + abs(err_t);
    select_indices(t) = find( sq_origin_dist == min(sq_origin_dist)  );
end

% for t = 1:T
%     b = B(:,t);
%     rel_err_t = err_select(1:M,t);
%     while rel_err_t(select_indices(t))/norm(B(:,t)) > 0.01
%         if select_indices(t) > 1
%             select_indices(t) = select_indices(t) - 1;
%         else
%             select_indices(t) = find(rel_err_t == min(rel_err_t));
%             break
%         end
%         selected_error(t) = err_select(select_indices(t),t);
%     end
% end



%% Plot AWMV
awmv_select = zeros(T,1);
for t = 1:T
    load(fullfile(output_dir,sprintf(baseFileName,select_indices(t),t)))
    load(fullfile(top_dir,om_name,dset_name,[P.prefix,'_',num2str(t),'.mat']))
    polar_vector = sum(polar_image,2);
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    az_signal = squeeze(sum(x_hat,1));
    var_sum = sum(az_signal(:));
    awmv_select(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end
figure(12)
plot(awmv)
title('AWMV Selected Indices')
save(fullfile(top_dir,[dset_name,om_name,'_mirror_indep_awmv.mat']),'awmv','select_indices','lambda_values','err_select','l1_select','rel_err_select','awmv')


    end
end