%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02';
mkdir(top_dir)

% Simulation name
sim_name = 'linearPoissonNoise';
file_name = 'lineSearchNoise';

% Define poisson dataset
N = 101;
T = 30;
[P,K,M,MM] = definePoissonP(N,T);
levels = 0.05:0.05:0.3;
alpha_vals = genPoissonScales(N,T,levels,'linear');

% Independent
% alg_name = 'IndepISM';
% indep_dir  = fullfile(top_dir,sim_name,alg_name);
% linearSimParamSearchPoisson(P,N,M,K,T,levels,alpha_vals,...
%                                 indep_dir,file_name)

% Coupled
alg_name = 'CoupledCGTV';
% coupled_dir  = fullfile(top_dir,sim_name,alg_name);       
% linearSimParamSearchCoupledPoisson(P,N,MM,K,T,levels,alpha_vals,...
%                                 coupled_dir,indep_dir,file_name)
  
%% Redo parameter selection

for nn = 1:numel(levels)
    awmv_all = zeros(MM,T);
    awmv_rmse = zeros(MM,1);
    mse_c = zeros(MM,1);
    l1_norm_c = zeros(MM,1);
    tv_penalty = zeros(MM,1);
    theta_stds1 = linspace(1,15,T);

    for i = 1:MM
        load(fullfile(top_dir,sim_name,alg_name,...
                [file_name,'_',num2str(nn),'.mat']))
        if i == 1
           A0 = unshifted_basis_vector_ft_stack_zpad(P);
        end
        for time = 1:T
            x = X_coupled(:,:,i,time);        
            fit = Ax_ft_1D(A0,x);
            mse_c(i) = mse_c(i) + norm( B(:,time)-fit ).^2/norm(B(:,time))^2;
            l1_norm_c(i) = P.params.lambda1(time)*l1_norm_c(i) + sum(abs(x(:)));
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_all(i,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        tv_penalty(i) = sum(abs(DiffPhiX_1D(X_coupled(:,:,i,:))),'all');
        awmv_rmse(i) = norm(awmv_all(i,:)-theta_stds1)/norm(theta_stds1);
    end

    crit1 = abs(mse_c-min(mse_c)).^2 /(max(mse_c)-min(mse_c))^2;
    crit2 = abs(l1_norm_c-min(l1_norm_c)).^2 /2/(max(l1_norm_c)-min(l1_norm_c))^2;
    crit3 = abs(tv_penalty-min(tv_penalty)).^2 /(max(tv_penalty)-min(tv_penalty))^2;
    crit = crit1+crit2+crit3;
    select_ind = find( (crit == min(crit)),1 );
    P.selected_lambda2 = P.lambda2_values(select_ind);
    P.coupled_select_ind = select_ind;
    save(fullfile(top_dir,sim_name,alg_name,...
                [file_name,'_',num2str(nn),'.mat']),'B','P','X_coupled')
end
                   
%% Plot results
awmv_coupled = zeros(numel(levels),T);
awmv_indep = zeros(numel(levels),T);
awmv_err_coupled = zeros(numel(levels),1);
awmv_err_indep = zeros(numel(levels),1);

% Coupled Error
alg_name = 'CoupledCGTV';
for nn = 1:numel(levels)
    load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'.mat']))
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
    end
    for time = 1:T
        x = X_coupled(:,:,P.coupled_select_ind,time);        
        fit_coupled = Ax_ft_1D(A0,x);
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_coupled(nn,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
    awmv_err_coupled(nn) = norm(theta_stds1-awmv_coupled(nn,:))/norm(theta_stds1);
end

% Independent Error
alg_name = 'IndepISM';
for nn = 1:numel(levels)
    load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'.mat']))
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
    end
    for time = 1:T
        x = X_indep(:,:,P.indep_select_ind(time),time);        
        fit_indep = Ax_ft_1D(A0,x);
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_indep(nn,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
    awmv_err_indep(nn) = norm(theta_stds1-awmv_indep(nn,:))/norm(theta_stds1);
end
     

figure(1)
hold on
plot(awmv_err_indep)
plot(awmv_err_coupled)
legend('indep','coupled')

figure(2)
for nn = 1:numel(levels)
    subplot(3,2,nn)
    hold on
    plot(theta_stds1)
    plot(awmv_indep(nn,:))
    plot(awmv_coupled(nn,:))
    legend('truth','indep','coupled','Location','Best')
end