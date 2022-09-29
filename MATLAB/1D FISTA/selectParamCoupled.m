function sel_ind = selectParamCoupled(MM,top_dir,sim_name,alg_name,ps_name,nn)
%UNTITLED10 Summary of this function goes here
mse_all = zeros(MM,1);
l1_norm_all = zeros(MM,1);
tv_penalty_all = zeros(MM,1);
awmv_coupled = zeros(T,numel(snr_levels),MM);
awmv_err_coupled = zeros(numel(snr_levels),MM);

for i = 1:MM
    load(fullfile(top_dir,sim_name,alg_name,...
        [ps_name,'_',num2str(nn),'_',num2str(i),'.mat']),...
        'mse','l1_norm','tv_penalty','P','X_hat','awmv_rmse')
    mse_all(i) = mse;
    l1_norm_all(i) = l1_norm;
    tv_penalty_all(i) = tv_penalty;
        awmv_coupled(:,nn,i) = computeAWMV_1D(X_hat,P.var_theta);
        awmv_err_coupled(nn,i) = norm(theta_stds-awmv_coupled(:,nn,i))/norm(theta_stds);
end
%     if nn == 1
%        A0 = unshifted_basis_vector_ft_stack_zpad(P);
%        theta_stds = P.theta_stds;
%     end
sel_ind = selectParamsCoupled(mse_all,l1_norm_all,tv_penalty_all);

for nn = 1:numel(snr_levels)
    [~,sel_ind(nn)] = min(awmv_err_coupled(nn,:));
end

end

