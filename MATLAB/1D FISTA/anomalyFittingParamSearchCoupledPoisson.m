%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
% top_dir = 'D:\CHESS_data\';
top_dir = '/cluster/home/dbanco02/data';


MM = 20;
noise_factor = 1;
NN = 1;

lambda2_values =  logspace(-0.3,2,MM);

num_ims = 30;
N = 101;
K = 20;
M = 50;
T = num_ims;
zPad = 0;
zMask = [];

theta_stds1 = [7*ones(1,T/2),12*ones(1,T/2)]';

%% Run ADMM
%{
for nn = 1
    % Input dirs
    dset_name = ['anomaly_noise_Poisson'];
    indep_name = '_indep_ISM';
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    indep_subdir = [dset_name,indep_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    indep_dir  = fullfile(top_dir,indep_subdir);

    mkdir(output_dir)
    load(fullfile(indep_dir,[dset_name,'_',num2str(M),'_','all']))

    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

    % Compute outputs
    rse_indep = zeros(M,T);
    mse_indep = zeros(M,T);
    l1_norm = zeros(M,T);
    for i = 1:M
        for time = 1:T
            x = X_indep(:,:,i,time);
            fit = Ax_ft_1D(A0ft_stack,x);
            l1_norm(i,time) = sum(abs(x(:)));
            mse_indep(i,time) = norm(fit-B(:,time));
            rse_indep(i,time) = norm(fit-B(:,time))/norm(B(:,time));   
        end
    end

    % Parameter Selection
    X_init = zeros(N,K,T);
    select_ind = zeros(T,1);
    lambda1_select = zeros(T,1);
    for time = 1:T
        crit = abs(l1_norm(:,time)*0.5).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
        lambda1_select(time) = P.lambda_values(select_ind(time));
        x = X_indep(:,:,select_ind(time),time);
        X_init(:,:,time) = x;
    end

    % Coupled Solution
    P.params.maxIter = 100;
    P.params.rho1 = 1.5;
    P.params.rho2 = 0.5;
    
    for i = 1:MM
        P.set = i;
        P.params.lambda1 = lambda1_select;
        P.params.lambda2 = lambda2_values(nn,i);
        [X_hat,~,~,~,~] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,B,X_init,P.params);
        save(fullfile(output_dir,[dset_name,'_',num2str(i),'_','time1']),...
            'B','X_hat','P');
    end
end

%}

%% Parameter Selection for coupled 

[ha_param, ~] = tight_subplot(4,3,[.005 .005],[.01 .01],[.01 .01]);
awmv_all = zeros(MM,T,NN);
awmv_rmse = zeros(MM,NN);
mse = zeros(MM,NN);
rse = zeros(MM,NN);
l1_norm = zeros(MM,NN);
tv_penalty = zeros(MM,NN);
for nn = 1:NN
    % Load coupled solution
    dset_name = ['anomaly_noise_Poisson'];
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    for i = 1:MM
        load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','time',num2str(1)]))
        if( (nn==1)&&(i==1) )
            % Construct dictionary
            A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
        end
        for time = 1:T
            x = X_hat(:,:,time);
            fit = Ax_ft_1D(A0ft_stack,x);
            mse(i,nn) = mse(i,nn) + norm( B(:,time)-fit ).^2;
            rse(i,nn) = rse(i,nn) + norm( B(:,time)-fit ).^2 /...
                                    norm( B(:,time))/T.^2;
            l1_norm(i,nn) = l1_norm(i,nn) + sum(abs(x(:)));
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_all(i,time,nn) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        tv_penalty(i,nn) = sum(abs(DiffPhiX_1D(X_hat)),'all');
        awmv_rmse(i,nn) = norm(awmv_all(i,:,nn)-theta_stds1')/norm(theta_stds1);
    end
    axes(ha_param(nn))
    plot(tv_penalty(:,nn), mse(:,nn),'-o')
end

figure(2)
plot(awmv_rmse)

%% Show L-curves and selected parameters
% % [~,s_i] = min(awmv_rmse);
% select_ind = zeros(NN,1);
% gamma_select = zeros(NN,T);
% figure(111)
% for nn = 1:NN
%     crit1 = abs(mse(:,nn)-min(mse(:,nn))).^2 /2/(max(mse(:,nn))-min(mse(:,nn)))^2;
%     crit2 = abs(l1_norm(:,nn)-min(l1_norm(:,nn))).^2 /2/(max(l1_norm(:,nn))-min(l1_norm(:,nn)))^2;
%     crit3 = abs(tv_penalty(:,nn)-min(tv_penalty(:,nn))).^2 /(max(tv_penalty(:,nn))-min(tv_penalty(:,nn)))^2;
%     crit = crit1+crit2+crit3;
%     select_ind(nn) = find( (crit == min(crit)),1 );
%     gamma_select(nn) = P.lambda_values(select_ind(nn));
%     % Plot L-curves
%     subplot(3,4,nn)
%     plot(crit1+crit2,...
%          crit3)
%     hold on
%     plot(crit1(select_ind(nn))+crit2(select_ind(nn)),...
%          crit3(select_ind(nn)),'or')
%     plot(crit1(s_i(nn))+crit2(s_i(nn)),...
%          crit3(s_i(nn)),'sb')
%     [awmv_rmse(select_ind(nn),nn),...
%      awmv_rmse(s_i(nn),nn)]
%     
%     hold off
%     xlabel('MSE')
%     ylabel('TV Penalty')
% 
% end

%% Show L-curves and selected parameters
% [~,s_i] = min(awmv_rmse);
select_ind = zeros(NN,1);
gamma_select = zeros(NN,T);
figure(111)
for nn = 1:NN
    crit1 = abs(mse(:,nn)-min(mse(:,nn))).^2 /2/(max(mse(:,nn))-min(mse(:,nn)))^2;
    crit2 = abs(l1_norm(:,nn)-min(l1_norm(:,nn))).^2 /2/(max(l1_norm(:,nn))-min(l1_norm(:,nn)))^2;
    crit3 = abs(tv_penalty(:,nn)-min(tv_penalty(:,nn))).^2 /(max(tv_penalty(:,nn))-min(tv_penalty(:,nn)))^2;
    crit = crit1+crit2+crit3;
    select_ind(nn) = find( (crit == min(crit)),1 );
    gamma_select(nn) = P.lambda_values(select_ind(nn));
    % Plot L-curves
    subplot(3,4,nn)
    plot(crit1+crit2,...
         crit3)
    hold on
    plot(crit1(select_ind(nn))+crit2(select_ind(nn)),...
         crit3(select_ind(nn)),'or')
%     plot(crit1(s_i(nn))+crit2(s_i(nn)),...
%          crit3(s_i(nn)),'sb')
%     [awmv_rmse(select_ind(nn),nn),...
%      awmv_rmse(s_i(nn),nn)]
    
    hold off
    xlabel('MSE')
    ylabel('TV Penalty')

end

%% Plot the pairs of AWMV fits 
% figure(222)
% for nn = 1:NN
%     subplot(3,4,nn)
%     hold on
%     plot(theta_stds1,'k')
%     plot(awmv_all(select_ind(nn),:,nn),'r')
%     plot(awmv_all(s_i(nn),:,nn),'b')
% end



%% 
% figure(22)
% hold on
% plot(theta_stds1,'Linewidth',2)
% for i = 1:2:MM
%     plot(awmv_all(i,:,11),'Linewidth',2)
% end


%% Show awmvs
close all
s_i2 = zeros(size(awmv_rmse));
f3 = figure(3);
f3.Position = [1000,500,600 400];
[ha_awmv, ~] = tight_subplot(4,3,[.005 .005],[.01 .01],[.01 .01]);
awmv_indep = zeros(T,NN);
awmv_rmse_indep = zeros(NN,1);
awmv_rmse_coupled = zeros(NN,1);
% awmv_rmse_coupled2 = zeros(NN,1);
for nn = 1:NN
%     awmv_rmse_coupled(nn) = awmv_rmse(s_i(nn),nn);
    % Load Independent and compute awmv
    dset_name = ['anomaly_noise_Poisson'];
    indep_name = '_indep_ISM';
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    indep_subdir = [dset_name,indep_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    indep_dir  = fullfile(top_dir,indep_subdir);
    load(fullfile(indep_dir,[dset_name,'_',num2str(M),'_','all']))
    if( (nn==1)&&(i==1) )
        % Construct dictionary
        A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
    end
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

    % Indep Parameter Selection
    X_init = zeros(N,K,T);
    select_ind = zeros(T,1);
    lambda1_select = zeros(T,1);
    for time = 1:T
        crit = abs(l1_norm(:,time)*0.50).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
        lambda1_select(time) = P.lambda_values(select_ind(time));
        x = X_indep(:,:,select_ind(time),time);
        X_init(:,:,time) = x;
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_indep(time,nn) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end
    
    % Coupled Parameter Selection
    crit = abs(tv_penalty(:,nn)).^2 + abs(mse(:,nn)).^2;
    s_i2(nn) = find( (crit == min(crit)),1 );
%     awmv_rmse_coupled2(nn) = awmv_rmse(s_i2(nn),nn);
    awmv_rmse_indep(nn) = norm(awmv_indep(:,nn)-theta_stds1)/norm(theta_stds1);
    axes(ha_awmv(nn))
    hold on
    plot(theta_stds1,'k--','Linewidth',2)
    plot(awmv_indep(:,nn),'b-','Linewidth',2)
    plot(awmv_all(s_i2(nn),:,nn),'r:','Linewidth',3) %s_i2(nn)
    NW = [min(xlim) max(ylim)]+[diff(xlim)*0.05 -diff(ylim)*0.1];
%     text(NW(1), NW(2), ['v=',sprintf('%1.2f',noise_std(nn))],'FontSize',14)
end
% Use Final plot for legend
axes(ha_awmv(NN+1))
hold on
plot(0,'k--','Linewidth',2)
plot(1,'b-','Linewidth',2)
plot(2,'r:','Linewidth',3)
legend('Truth','\gamma = 0','\gamma = \gamma*',...
       'FontSize',16,'EdgeColor',[1 1 1],'location','Northwest')


% PLot RMSE in AWMV
figure(4)
hold on
plot(awmv_rmse_indep,'b','Linewidth',2)
plot(awmv_rmse_coupled,'r','Linewidth',2)
% plot(awmv_rmse_coupled2,'g')
legend('\gamma = 0','\gamma = \gamma*','Location','Northwest','FontSize',16)
ylabel('Relative MSE in AWMV','FontSize',16)
xlabel('\sigma_{noise}','FontSize',16)


