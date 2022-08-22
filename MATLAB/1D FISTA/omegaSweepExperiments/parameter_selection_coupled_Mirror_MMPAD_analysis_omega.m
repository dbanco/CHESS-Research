clear all
close all

% Parameter selection
disp('Setup params')
P.set = 1;

% Parent directory
% top_dir = 'E:\MMPAD_omega';
top_dir = '/cluster/home/dbanco02/data/MMPAD_omega';
om_dir = {'omega2','omega3','omega4','omega5'};
r_dir = {'ring1','ring2','ring3','ring4'};

for o = 1:3
for ring_num = 1:4
% Input dirs
dset_name = r_dir{ring_num};
om_name = om_dir{o};

% Output dirs
output_name = '_coupled_CG_TVphi_Mirror';
output_subdir = [dset_name,om_dir{o},output_name];

% Setup directories
dataset =  fullfile(top_dir,om_dir{o},dset_name);
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
    e_data = load(fullfile(output_dir,sprintf(baseFileName,i)),'P','X_hat');
    
    tv_penalty(i) = sum(abs(DiffPhiX_1D(e_data.X_hat)),'all');
    
    for j = 1:T
        b = B(:,j);
        x = squeeze(e_data.X_hat(:,:,j));
        fit = Ax_ft_1D(A0ft_stack,x);   
        err_select(i,j) = sum(( fit(:) - b(:) ).^2)/norm(b)^2;
        l0_select(i,j) = sum(x(:) > 1e-4*sum(x(:)));
%         l1_select(i,j) = sum(x(:));
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        vdfs(:,i,j) = az_signal./var_sum;
        awmv_az(i,j) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        
        if j == 40
            hold off
            figure(1)
            plot(b)
            hold on
            plot(fit)
            legend('data','fit')
        end
    end
    im_ind = im_ind + 1;
end

%% Criterion separate params
% noise_eta = 0.2;
% discrep_crit = abs(err_select'-noise_eta);
% 
% [lambda_indices,~] = find(discrep_crit' == min(discrep_crit'));
% param_select = P.lambda_values(lambda_indices);
% lambda_values_separate = param_select;
% 
% % Criterion single param
% discrep_crit = abs(mean(err_select,2)-noise_eta);
% lambda_index = find(discrep_crit == min(discrep_crit));
% param_select_single = P.lambda_values(lambda_index);
% lambda_values_single = ones(T,1)*param_select_single;

% figure(2)
% plot(param_select,'o-')
% title('Parameters selected')


%% Plot
lambda_vals = P.lambda_values;

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
% sq_origin_dist = total_l1.^2 + total_err.^2;
% select_ind = find(sq_origin_dist == min(sq_origin_dist));
err_scale = total_err/max(total_err(:));
tv_scale = tv_penalty/max(tv_penalty(:));
sq_origin_dist = abs(tv_scale).^2 + abs(err_scale).^2;
select_ind = find(sq_origin_dist == min(sq_origin_dist),1);

save([dset_name,om_name,'_mirror_coupled_awmv.mat'],'awmv_az','select_ind',...
    'lambda2_values','err_select','l1_select','tv_penalty')
end
end
%{
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
%}

