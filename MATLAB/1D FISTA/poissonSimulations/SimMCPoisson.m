function SimMCPoisson(P,N,K,T,levels,alpha_vals,...
                 indepPS_dir,ps_name,mc_dir,sim)

% Fits for different parameters/noise levels
for nn = 1:numel(levels)
    
    % Make folder to save trials
    trials_dir = fullfile(indepPS_dir,[mc_dir,'_',num2str(nn)]);
    mkdir(trials_dir)
    
    % Get selected parameter indices
    ps = load(fullfile(indepPS_dir,[ps_name,'_',num2str(nn),'.mat']),'P');
    lambdas = ps.P.selected_lambdas;
       
    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

    % Independent Solution
    x_init = zeros(N,K);

    for i = 1:P.trials
        X_indep = zeros(N,K,T);
        awmv = zeros(T,1);
        [B,~,~,~] = genSimDataPoisson(N,T,alpha_vals(nn),sim);
        P.set = i;
        for t = 1:T
            % Solve
             P.params.lambda1 = lambdas(t);
            [x_hat] = convADMM_LASSO_Sherman_1D(A0ft_stack,B(:,t),x_init,P.params);  
            X_indep(:,:,t) = x_hat;
            az_signal = squeeze(sum(x_hat,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        save(fullfile(trials_dir,['trial_',num2str(i)]),...
             'B','X_indep','P','awmv');
    end
end
