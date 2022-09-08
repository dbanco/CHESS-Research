function SimMCPoissonCoupled(P,N,K,T,levels,alpha_vals,...
                     indepPS_dir,coupledPS_dir,ps_name,mc_dir,sim)

                    % Fits for different parameters/noise levels
NN = numel(levels);
for nn =  6%1:NN
    
    % Make folder to save trials
    indepTrials_dir = fullfile(indepPS_dir,[mc_dir,'_',num2str(nn)]);
    coupledTrials_dir = fullfile(coupledPS_dir,[mc_dir,'_',num2str(nn)]);
    mkdir(coupledTrials_dir)
    
    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
    
    % Load coupled selected parameter (selected_lambda2)
    ps1 = load( fullfile(coupledPS_dir,[ps_name,'_',num2str(nn)]) );
    ps2 = load(fullfile(indepPS_dir,[ps_name,'_',num2str(nn),'.mat']),'P');
    
	for i = 1:P.trials
        % Load indep solution and data
        indepData = load(fullfile(indepTrials_dir,['trial_',num2str(i)]));
        B = indepData.B;
        
        % Coupled Solution
        P.params.maxIter = 100;
        P.params.rho1 = 1.5;
        P.params.rho2 = 0.5;
        P.set = i;
        P.params.lambda1 = ps2.P.selected_lambdas;
        P.params.lambda2 = ps1.P.selected_lambda2;
        % Solve
        [X_coupled,~,~,~,~] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,indepData.B,...
                                               indepData.X_indep,P.params);                              

        awmv = zeros(T,1);
        for time = 1:T
            x = X_coupled(:,:,time);        
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        save(fullfile(coupledTrials_dir,['trial_',num2str(i)]),...
            'B','X_coupled','P','awmv');
	end
end
end
