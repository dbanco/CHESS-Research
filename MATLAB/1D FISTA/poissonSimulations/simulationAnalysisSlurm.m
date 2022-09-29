top_dir = 'D:\Simulations';
sim = 'linear';
sim_name = [sim,'PoissonNoise6'];
file_name = 'lineSearchNoise';
N = 101;
T = 30;
[P,K,M,MM] = definePoissonP(N,T,sim);
theta_stds = P.theta_stds;
snr_levels = [10,5,2.5,2,1.75,1.5,1.25,1.1];

close all

awmv_coupled = zeros(T,numel(snr_levels),MM);
awmv_indep = zeros(T,numel(snr_levels));
awmv_err_coupled = zeros(numel(snr_levels),MM);
awmv_err_indep = zeros(numel(snr_levels),1);
sel_ind = zeros(numel(snr_levels),1);
awmv_err_coupled_sel = zeros(numel(snr_levels),1);
% Coupled Error
alg_name = 'CoupledCGTV';
for nn = 1:numel(snr_levels)
    mse_all = zeros(MM,1);
    l1_norm_all = zeros(MM,1);
    tv_penalty_all = zeros(MM,1);
    for i = 1:MM
        load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'_',num2str(i),'.mat']),...
            'mse','l1_norm','tv_penalty','P','X_hat','awmv_rmse')
        mse_all(i) = mse;
        l1_norm_all(i) = l1_norm;
        tv_penalty_all(i) = tv_penalty;
        awmv_coupled(:,nn,i) = computeAWMV_1D(X_hat,P.var_theta);
        awmv_err_coupled(nn,i) = norm(theta_stds-awmv_coupled(:,nn,i))/norm(theta_stds);
    end
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
       theta_stds = P.theta_stds;
    end
    sel_ind(nn) = selectParamsCoupled(mse_all,l1_norm_all,tv_penalty_all);
end
for nn = 1:numel(snr_levels)
    [awmv_err_coupled_sel(nn),sel_ind(nn)] = min(awmv_err_coupled(nn,:));
end

% Independent Error
alg_name = 'IndepISM';
for nn = 1:numel(snr_levels)
    load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'.mat']))
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
    end
    for time = 1:T
        x = X_indep(:,:,P.indep_select_ind(time),time);        
        awmv_indep(time,nn) = computeAWMV_1D(x,P.var_theta);
    end
    
    awmv_err_indep(nn) = norm(theta_stds-awmv_indep(:,nn))/norm(theta_stds);
end

figure(1)
hold on
plot(awmv_err_indep)
% plot(awmv_err_coupled(:,30))
plot(awmv_err_coupled_sel)
legend('indep','coupled')

figure(2)
for nn = 1:numel(snr_levels)
    subplot(3,3,nn)
    hold on
    plot(theta_stds)
    plot(awmv_indep(:,nn))
    plot(awmv_coupled(:,nn,sel_ind(nn)))
    legend('truth','indep','coupled','Location','Best')
end


