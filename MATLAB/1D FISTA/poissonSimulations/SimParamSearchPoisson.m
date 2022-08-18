function SimParamSearchPoisson(P,N,M,K,T,levels,...
                            alpha_vals,output_dir,file_name,sim)
mkdir(output_dir)

% Fits for different parameters/noise levels
for nn = 1:numel(levels)
    % Setup directories
    f_name =  [file_name,'_',num2str(nn),'.mat'];
       
    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

    % Generate data
    if strcmp(sim,'linear')
        [B,~,theta_stds1,~] = genLinearPoisson(N,T,alpha_vals(nn));
    elseif strcmp(sim,'anomaly')
        [B,~,theta_stds1,~] = genAnomalyPoisson(N,T,alpha_vals(nn));
    end

    % Independent Solution
    x_init = zeros(N,K);
    X_indep = zeros(N,K,M,T);
    for i = 1:M
        P.set = i;
        P.params.lambda1 = P.lambda_values(i);
        for t = 1:T
            % Solve
            [x_hat] = convADMM_LASSO_Sherman_1D(A0ft_stack,B(:,t),x_init,P.params);  
            X_indep(:,:,i,t) = x_hat;
        end
    end
    
    % Parameter selection
    [mse_indep,l1_norm,~,~,~] = exploreParametersIndep(X_indep,P,B);
    select_ind = selectParamsIndep(mse_indep,l1_norm);
    
    P.indep_select_ind = select_ind;
    P.selected_lambdas = P.lambda_values(select_ind);
    
    save(fullfile(output_dir,f_name),'B','X_indep','P');
end
