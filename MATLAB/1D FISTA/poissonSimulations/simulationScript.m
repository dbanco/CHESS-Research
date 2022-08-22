%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/home/dbanco02/data/Simulations';
mkdir(top_dir)

% Simulation name
sim = 'anomaly';
sim_name = [sim,'PoissonNoise2'];
file_name = 'lineSearchNoise';

% Define poisson dataset
N = 101;
T = 30;
[P,K,M,MM] = definePoissonP(N,T);
levels = 0.05:0.05:0.3;
[alpha_vals,theta_stds] = genPoissonScales(N,T,levels,sim);
P.theta_stds = theta_stds;

%% Run algorithm

% Independent
alg_name = 'IndepISM';
indep_dir  = fullfile(top_dir,sim_name,alg_name);
% SimParamSearchPoisson(P,N,M,K,T,levels,alpha_vals,...
%                                 indep_dir,file_name,sim)

% Coupled
alg_name = 'CoupledCGTV';
coupled_dir  = fullfile(top_dir,sim_name,alg_name);       
% SimParamSearchCoupledPoisson(P,N,MM,K,T,levels,alpha_vals,...
%                                 coupled_dir,indep_dir,file_name,sim)

%% Redo parameter selection indep

% alg_name = 'IndepISM';
% indep_dir = fullfile(top_dir,sim_name,alg_name);
% for nn = 1:numel(levels)
%     f_name = [file_name,'_',num2str(nn),'.mat'];
%     load(fullfile(indep_dir,f_name));
% 
%     [mse_indep,l1_norm,awmv,~,~] = exploreParametersIndep(X_indep,P,B);
%     select_ind = selectParamsIndep(mse_indep,l1_norm);
% %     select_ind = selectParamsIndepAWMV(awmv,theta_stds);
%     awmv_select = selectAWMV(awmv,select_ind);
%     P.indep_select_ind = select_ind;
%     P.selected_lambdas = P.lambda_values(select_ind);
% %     save(fullfile(indep_dir,f_name),'B','X_indep','P');
%     fprintf('%i, ',nn)
% end
% close all
% figure(3)
% hold on
% plot(awmv_select)
% plot(theta_stds)

%% Redo parameter selection coupled
alg_name = 'CoupledCGTV';
coupled_dir  = fullfile(top_dir,sim_name,alg_name); 
for nn = 1:numel(levels)
    load(fullfile(coupled_dir,[file_name,'_',num2str(nn),'.mat']))
    [mse,l1_norm,tv_penalty,awmv,~,B] = exploreParametersCoupled(X_coupled,P,B);
%     select_ind = selectParamsCoupled(mse,l1_norm,tv_penalty);
    select_ind = selectParamsCoupledAWMV(awmv,theta_stds);

%     awmv_rmse = zeros(MM,1);
%     awmv_rmse(i) = norm(awmv(i,:)-theta_stds1)/norm(theta_stds1);

    P.selected_lambda2 = P.lambda2_values(select_ind);
    P.coupled_select_ind = select_ind;
    save(fullfile(coupled_dir,...
                [file_name,'_',num2str(nn),'.mat']),'B','P','X_coupled')
            fprintf('%i, ',nn)
end  

close all
figure(3)
hold on
plot(awmv(:,select_ind))
plot(theta_stds)
%% Plot results
% awmv_coupled = zeros(numel(levels),T);
% awmv_indep = zeros(numel(levels),T);
% awmv_err_coupled = zeros(numel(levels),1);
% awmv_err_indep = zeros(numel(levels),1);
% hand_ind = [15,17,20,38,35,50];
% % hand_ind = [15,17,20,5,40,40];
% 
% % Coupled Error
% alg_name = 'CoupledCGTV';
% for nn = 1:numel(levels)
%     load(fullfile(top_dir,sim_name,alg_name,...
%             [file_name,'_',num2str(nn),'.mat']))
%     if nn == 1
%        A0 = unshifted_basis_vector_ft_stack_zpad(P);
% %        theta_stds = P.theta_stds;
%         theta_stds = theta_stds;
%     end
%     for time = 1:T
%         x = X_coupled(:,:,hand_ind(nn),time);        
%         fit_coupled = Ax_ft_1D(A0,x);
%         az_signal = squeeze(sum(x,1));
%         var_sum = squeeze(sum(az_signal(:)));
%         awmv_coupled(nn,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     end
%     awmv_err_coupled(nn) = norm(theta_stds-awmv_coupled(nn,:))/norm(theta_stds);
% end
% 
% % Independent Error
% alg_name = 'IndepISM';
% for nn = 1:numel(levels)
%     load(fullfile(top_dir,sim_name,alg_name,...
%             [file_name,'_',num2str(nn),'.mat']))
%     if nn == 1
%        A0 = unshifted_basis_vector_ft_stack_zpad(P);
%     end
%     for time = 1:T
%         x = X_indep(:,:,P.indep_select_ind(time),time);        
%         fit_indep = Ax_ft_1D(A0,x);
%         az_signal = squeeze(sum(x,1));
%         var_sum = squeeze(sum(az_signal(:)));
%         awmv_indep(nn,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
%     end
%     awmv_err_indep(nn) = norm(theta_stds-awmv_indep(nn,:))/norm(theta_stds);
% end
%      
% 
% figure(1)
% hold on
% plot(awmv_err_indep)
% plot(awmv_err_coupled)
% legend('indep','coupled')
% 
% figure(2)
% for nn = 1:numel(levels)
%     subplot(3,2,nn)
%     hold on
%     plot(theta_stds)
%     plot(awmv_indep(nn,:))
%     plot(awmv_coupled(nn,:))
%     legend('truth','indep','coupled','Location','Best')
% end

%% Plot results
awmv_coupled = zeros(numel(levels),T);
awmv_indep = zeros(numel(levels),T);
awmv_err_coupled = zeros(numel(levels),1);
awmv_err_indep = zeros(numel(levels),1);
% Coupled Error
alg_name = 'CoupledCGTV';
for nn = 1:numel(levels)
    load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'.mat']))
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
%        theta_stds = P.theta_stds;
        theta_stds = theta_stds;
    end
    for time = 1:T
        x = X_coupled(:,:,P.coupled_select_ind,time);        
        fit_coupled = Ax_ft_1D(A0,x);
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_coupled(nn,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
    awmv_err_coupled(nn) = norm(theta_stds-awmv_coupled(nn,:))/norm(theta_stds);
    P.coupled_select_ind
end

% Independent Error
alg_name = 'IndepISM';
for nn = 1:numel(levels)
    load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'.mat']))
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
    end
    for time = 1:T
        x = X_indep(:,:,P.indep_select_ind(time),time);        
        fit_indep = Ax_ft_1D(A0,x);
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_indep(nn,time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
    awmv_err_indep(nn) = norm(theta_stds-awmv_indep(nn,:))/norm(theta_stds);
end
     

figure(1)
hold on
plot(awmv_err_indep)
plot(awmv_err_coupled)
legend('indep','coupled')

figure(2)
for nn = 1:numel(levels)
    subplot(3,2,nn)
    hold on
    plot(theta_stds)
    plot(awmv_indep(nn,:))
    plot(awmv_coupled(nn,:))
    legend('truth','indep','coupled','Location','Best')
end

%% Plot awmv for different parameter values
awmv_coupled = zeros(numel(levels),T,MM);
awmv_indep = zeros(numel(levels),T,1);
awmv_err_coupled = zeros(numel(levels),MM);
awmv_err_indep = zeros(numel(levels),1);

% Coupled Error
alg_name = 'CoupledCGTV';
for nn = 1:numel(levels)
    load(fullfile(top_dir,sim_name,alg_name,...
            [file_name,'_',num2str(nn),'.mat']))
    if nn == 1
       A0 = unshifted_basis_vector_ft_stack_zpad(P);
    end
    for time = 1:T
        for ii = 1:MM
            x = X_coupled(:,:,ii,time);        
            fit_coupled = Ax_ft_1D(A0,x);
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_coupled(nn,time,ii) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
    end
    awmv_err_coupled(nn,ii) = norm(theta_stds-awmv_coupled(nn,:,ii))/norm(theta_stds);
end

%% Explore indep parameters
alg_name = 'IndepISM';
indep_dir = fullfile(top_dir,sim_name,alg_name);
iii = [30,20,20,18,17,17];
for nn = 1:6
f_name = [file_name,'_',num2str(nn),'.mat'];
load(fullfile(indep_dir,f_name));
[mse_indep,l1_norm,awmv_indep,fits,B] = exploreParametersIndep(X_indep,P,B);

ii = iii(nn);
% figure(1)
% hold on
% plot(theta_stds1)
% plot(awmv_indep(:,ii))

figure(2)
subplot(6,1,nn)
t = 30;
hold on
plot(B(:,t),'Linewidth',2)
plot(fits(:,t,25))
legend('data','recon')
end