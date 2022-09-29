function sel_ind = selectParamCoupled(MM,theta_stds,top_dir,sim_name,alg_name,ps_name,nn)
%UNTITLED10 Summary of this function goes here
mse_all = zeros(MM,1);
l1_norm_all = zeros(MM,1);
tv_penalty_all = zeros(MM,1);
awmv_err_coupled = zeros(MM,1);

for i = 1:MM
    load(fullfile(top_dir,sim_name,alg_name,...
        [ps_name,'_',num2str(nn),'_',num2str(i),'.mat']),...
        'mse','l1_norm','tv_penalty','P','X_hat','X_hat')
    mse_all(i) = mse;
    l1_norm_all(i) = l1_norm;
    tv_penalty_all(i) = tv_penalty;
        awmv = computeAWMV_1D(X_hat,P.var_theta);
        awmv_err_coupled(i) = norm(theta_stds(:)-awmv(:))/norm(theta_stds(:));
end
sel_ind = selectParamsCoupled(mse_all,l1_norm_all,tv_penalty_all);
[~,sel_ind] = min(awmv_err_coupled);


end

