function SimParamSearchCoupledPoisson2(P,MM,snr_levels,...
                  output_dir,indep_dir,file_name,sim_name)
mkdir(output_dir)

% Fits for different parameters/noise levels
N = P.num_theta;
K = P.num_var_t;
T = P.num_ims;
NN = numel(snr_levels);
for nn = 8:NN
    % Setup directories
    f_name =  [file_name,'_',num2str(nn),'.mat'];
       
    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
    
    % Load indep solution
    i_name = [file_name,'_',num2str(nn),'.mat'];
    indepData = load(fullfile(indep_dir,i_name));
    P.params.lambda1 = indepData.P.selected_lambdas;
    x_init = zeros(N,K,T);
    for t = 1:T
        ii = indepData.P.indep_select_ind(t);
        x_init(:,:,t) = squeeze(indepData.X_indep(:,:,ii,t));
    end
    % Coupled Solution
    P.params.maxIter = 100;
    P.params.rho1 = 1;
    P.params.rho2 = 1;
    X_coupled = zeros(N,K,MM,T);
    for i = 47:MM 
        P.set = i;
        P.params.lambda2 = P.lambda2_values(i);
        % Solve
        [X_hat,~,~,~,~] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,indepData.B,...
                                                     x_init,P.params); 
        X_coupled(:,:,i,:) = X_hat;                             
    end
    
    awmv_all = zeros(MM,T);
    awmv_rmse = zeros(MM,1);
    mse_c = zeros(MM,1);
    l1_norm_c = zeros(MM,1);
    tv_penalty = zeros(MM,1);
    for i = 1:MM
        for time = 1:T
            x = X_coupled(:,:,i,time);        
            fit = Ax_ft_1D(A0ft_stack,x);
            mse_c(i) = mse_c(i) + norm( indepData.B(:,time)-fit ).^2/norm(indepData.B(:,time))^2;
            l1_norm_c(i) = P.params.lambda1(time)*l1_norm_c(i) + sum(abs(x(:)));
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_all(i,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        tv_penalty(i) = sum(abs(DiffPhiX_1D(squeeze(X_coupled(:,:,i,:)))),'all');
        awmv_rmse(i) = norm(awmv_all(i,:)-P.theta_stds)/norm(P.theta_stds);
    end
    
    crit1 = abs(mse_c-min(mse_c)).^2 /(max(mse_c)-min(mse_c))^2;
    crit2 = abs(l1_norm_c-min(l1_norm_c)).^2 /2/(max(l1_norm_c)-min(l1_norm_c))^2;
    crit3 = abs(tv_penalty-min(tv_penalty)).^2 /(max(tv_penalty)-min(tv_penalty))^2;
    crit = crit1+crit2+crit3;
    select_ind = find( (crit == min(crit)),1 );
    
    P.coupled_select_ind = select_ind;
    P.selected_lambda2 = P.lambda2_values(select_ind);
    
    save(fullfile(output_dir,f_name),'B','X_coupled','P');
end
