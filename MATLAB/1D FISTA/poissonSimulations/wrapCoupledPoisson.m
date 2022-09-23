function wrapCoupledPoisson(P)

N = P.num_theta;
T = P.num_ims;
K = P.num_var_t;
i = P.set;
P.params.lambda2 = P.lambda2_values(i);

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
            
% Load indep solution to get parameters
indepData = load(P.indepFile);
P.params.lambda1 = indepData.P.selected_lambdas;
x_init = zeros(N,K,T);
for t = 1:T
    ii = indepData.P.indep_select_ind(t);
    x_init(:,:,t) = squeeze(indepData.X_indep(:,:,ii,t));
end

% Solve
[X_hat,~,~,~,~] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,indepData.B,...
                                             x_init,P.params);      

% Compute outputs                              
awmv = computeAWMV_1D(X_hat,P.var_theta);
mse = 0; l1_norm = 0;
for time = 1:T
    x = X_hat(:,:,time);        
    fit = Ax_ft_1D(A0ft_stack,x);
    mse = mse + norm( indepData.B(:,time)-fit ).^2;
    l1_norm = l1_norm + P.params.lambda1(time)*sum(abs(x(:)));
end
tv_penalty = sum(abs(DiffPhiX_1D(X_hat)),'all');
awmv_rmse = norm(awmv-P.theta_stds)/norm(P.theta_stds);

save(P.outputFile,'X_hat','P','mse','l1_norm','tv_penalty','awmv_rmse');

end

