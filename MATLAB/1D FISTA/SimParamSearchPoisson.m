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
    
    % Error/L1-norm
    mse_indep = zeros(M,T);
    l1_norm = zeros(M,T);
    for i = 1:M
        for time = 1:T
            x = X_indep(:,:,i,time);
            fit = Ax_ft_1D(A0ft_stack,x);
            l1_norm(i,time) = sum(abs(x(:)));
            mse_indep(i,time) = norm(fit-B(:,time));
        end
    end
    
    % Parameter selection
    select_ind = zeros(T,1);
    for time = 1:T
        crit = abs(l1_norm(:,time)*0.5).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
    end
    
    P.indep_select_ind = select_ind;
    P.selected_lambdas = P.lambda_values(select_ind);
    
    save(fullfile(output_dir,f_name),'B','X_indep','P');
end
