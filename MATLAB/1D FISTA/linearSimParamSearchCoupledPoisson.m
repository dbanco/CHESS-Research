function linearSimParamSearchCoupledPoisson(P,N,MM,K,T,levels,...
                  alpha_vals,output_dir,indep_dir,file_name)
mkdir(output_dir)

% Fits for different parameters/noise levels
NN = numel(levels);
for nn = 1:NN
    % Setup directories
    f_name =  [file_name,'_',num2str(nn),'.mat'];
       
    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

    % Generate data
    [B,theta_stds1] = genLinearPoisson(N,T,alpha_vals(nn));

    % Load indep solution
    i_name = [file_name,'_',num2str(nn),'.mat'];
    indepData = load(fullfile(indep_dir,i_name));

    
    % Coupled Solution
    P.params.maxIter = 100;
    P.params.rho1 = 1.5;
    P.params.rho2 = 0.5;
    X_coupled = zeros(N,K,MM,T);
    for i = 1:MM 
        P.set = i;
        P.params.lambda1 = indepData.P.selected_lambdas;
        P.params.lambda2 = P.lambda2_values(i);
        % Solve
        x_init = squeeze(indepData.X_indep(:,:,i,:));
        [X_hat,~,~,~,~] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,B,...
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
            x = X_hat(:,:,time);        
            fit = Ax_ft_1D(A0ft_stack,x);
            mse_c(i) = mse_c(i) + norm( B(:,time)-fit ).^2/norm(B(:,time))^2;
            l1_norm_c(i) = P.params.lambda1(time)*l1_norm_c(i) + sum(abs(x(:)));
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_all(i,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        tv_penalty(i) = sum(abs(DiffPhiX_1D(X_hat)),'all');
        awmv_rmse(i) = norm(awmv_all(i,:)-theta_stds1')/norm(theta_stds1);
    end
    
    crit1 = abs(mse_c(:,nn)-min(mse_c(:,nn))).^2 /2/(max(mse_c(:,nn))-min(mse_c(:,nn)))^2;
    crit2 = abs(l1_norm_c(:,nn)-min(l1_norm_c(:,nn))).^2 /2/(max(l1_norm_c(:,nn))-min(l1_norm_c(:,nn)))^2;
    crit3 = abs(tv_penalty(:,nn)-min(tv_penalty(:,nn))).^2 /(max(tv_penalty(:,nn))-min(tv_penalty(:,nn)))^2;
    crit = crit1+crit2+crit3;
    select_ind = find( (crit == min(crit)),1 );
    
    P.coupled_select_ind = select_ind;
    P.selected_lambda2 = P.lambda2_values(select_ind);
    
    save(fullfile(output_dir,f_name),'B','X_coupled','P');
end
