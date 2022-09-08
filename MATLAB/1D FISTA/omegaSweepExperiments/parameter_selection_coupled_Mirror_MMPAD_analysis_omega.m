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

for o = 1:4
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

load([dset_name,om_name,'_mirror_coupled_awmv.mat'])
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

