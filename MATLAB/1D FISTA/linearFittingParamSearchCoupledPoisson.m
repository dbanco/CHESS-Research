%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\CHESS_data\';
% top_dir = '/cluster/shared/dbanco02';

noise_factor = [0.1:0.1:1];
NN = numel(noise_factor);

MM = 15;

lambda2_values = [logspace(-2.2,-1,MM);
                  logspace(-2,-0.8,MM);
                  logspace(-2,-0.8,MM);
                  logspace(-2,-0.5,MM);
                  logspace(-1,-0.5,MM);
                  logspace(-1,-0.5,MM);
                  logspace(-1,-0.5,MM);
                  logspace(-1,-0.2,MM);
                  logspace(-1,-0.2,MM);
                  logspace(-1,-0.2,MM);
                  logspace(-1,-0.2,MM)];

num_ims = 50;
N = 101;
K = 20;
M = 50;
T = num_ims;
zPad = 0;
zMask = [];

theta_stds1 = linspace(1,15,T)';

%% Run ADMM

for nn = 1:NN
    % Input dirs
    dset_name = ['singlePeak_noise_poisson',num2str(nn)];
    indep_name = '_indep_ISM';
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    indep_subdir = [dset_name,indep_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    indep_dir  = fullfile(top_dir,indep_subdir);

    mkdir(output_dir)
    load(fullfile(indep_dir,[dset_name,'_',num2str(num_ims),'_','all']))

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


%% Parameter Selection for coupled 
[ha_param, ~] = tight_subplot(4,3,[.005 .005],[.01 .01],[.01 .01]);
awmv_all = zeros(MM,T,NN);
awmv_rmse = zeros(MM,NN);
mse_c = zeros(MM,NN);
l1_norm_c = zeros(MM,NN);
tv_penalty = zeros(MM,NN);
for nn = 1:NN
    % Load coupled solution
    dset_name = ['singlePeak_noise_poisson',num2str(nn)];
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
            mse_c(i,nn) = mse_c(i,nn) + norm( B(:,time)-fit ).^2/norm(B(:,time))^2;
%             figure(898)
%             plot(B(:,time),'Linewidth',2)
%             hold on
%             plot(fit)
%             title(sprintf('RSE = %0.3f',norm( B(:,time)-fit ).^2/norm(B(:,time))^2))
%             hold off
%             pause()
            l1_norm_c(i,nn) = P.params.lambda1(time)*l1_norm_c(i,nn) + sum(abs(x(:)));
            az_signal = squeeze(sum(x,1));
            var_sum = squeeze(sum(az_signal(:)));
            awmv_all(i,time,nn) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
        end
        tv_penalty(i,nn) = sum(abs(DiffPhiX_1D(X_hat)),'all');
        awmv_rmse(i,nn) = norm(awmv_all(i,:,nn)-theta_stds1')/norm(theta_stds1);
    end
    axes(ha_param(nn))
    plot(tv_penalty(:,nn), mse_c(:,nn),'-o')
end

figure(2)
plot(awmv_rmse)
%% Show awmvs
close all
[~,s_i] = min(awmv_rmse);
s_i2 = zeros(NN,1);
gamma_select = zeros(NN,T);
for nn = 1:NN
    crit1 = abs(mse_c(:,nn)-min(mse_c(:,nn))).^2 /2/(max(mse_c(:,nn))-min(mse_c(:,nn)))^2;
    crit2 = abs(l1_norm_c(:,nn)-min(l1_norm_c(:,nn))).^2 /2/(max(l1_norm_c(:,nn))-min(l1_norm_c(:,nn)))^2;
    crit3 = abs(tv_penalty(:,nn)-min(tv_penalty(:,nn))).^2 /(max(tv_penalty(:,nn))-min(tv_penalty(:,nn)))^2;
    crit = crit1+crit2+crit3;
    s_i2(nn) = find( (crit == min(crit)),1 );
end
f3 = figure(3);
f3.Position = [1000,500,600 400];
[ha_awmv, ~] = tight_subplot(4,3,[.005 .005],[.01 .01],[.01 .01]);
awmv_indep = zeros(T,NN);
awmv_rmse_indep = zeros(NN,1);
awmv_rmse_coupled = zeros(NN,1);
for nn = 1:NN
    rms = sqrt(sum(B(:,1).^2)/N);
%     snr(nn) = rms/noise_std(nn);
    awmv_rmse_coupled(nn) = awmv_rmse(s_i2(nn),nn);
    % Load Independent and compute awmv
    dset_name = ['singlePeak_noise',num2str(nn)];
    indep_name = '_indep_ISM';
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    indep_subdir = [dset_name,indep_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    indep_dir  = fullfile(top_dir,indep_subdir);
    load(fullfile(indep_dir,[dset_name,'_',num2str(num_ims),'_','all2']))
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

    % Parameter Selection
    X_init = zeros(N,K,T);
    select_ind = zeros(T,1);
    lambda1_select = zeros(T,1);
    for time = 1:T
        crit = abs(l1_norm(:,time)*0.45).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
        lambda1_select(time) = P.lambda_values(select_ind(time));
        x = X_indep(:,:,select_ind(time),time);
        X_init(:,:,time) = x;
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_indep(time,nn) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;        
    end
    
    awmv_rmse_indep(nn) = norm(awmv_indep(:,nn)-theta_stds1)/norm(theta_stds1);
    awmv_noise = awmv_indep(:,nn)-theta_stds1
    figure(88)
    hist(awmv_noise-mean(awmv_noise))

    axes(ha_awmv(nn))
    hold on
    plot(theta_stds1,'k--','Linewidth',2)
    plot(awmv_indep(:,nn),'b-','Linewidth',2)
    plot(awmv_all(s_i2(nn),:,nn),'r:','Linewidth',3)
    NW = [min(xlim) max(ylim)]+[diff(xlim)*0.05 -diff(ylim)*0.1];
    text(NW(1), NW(2), ['v=',sprintf('%1.2f',noise_factor(nn))],'FontSize',14)
end
% Use Final plot for legend
axes(ha_awmv(NN+1))
hold on
plot(0,'k--','Linewidth',2)
plot(1,'b-','Linewidth',2)
plot(2,'r:','Linewidth',3)
legend('Truth','\gamma = 0','\gamma = \gamma*',...
       'FontSize',16,'EdgeColor',[1 1 1],'location','Northwest')

   %% Show L-curves and selected parameters
nn = 3;

[~,s_i] = min(awmv_rmse);
select_ind = zeros(NN,1);
gamma_select = zeros(NN,T);
figure(111)
for nn = 1:NN
    crit1 = abs(mse_c(:,nn)-min(mse_c(:,nn))).^2 /2/(max(mse_c(:,nn))-min(mse_c(:,nn)))^2;
    crit2 = abs(l1_norm_c(:,nn)-min(l1_norm_c(:,nn))).^2 /2/(max(l1_norm_c(:,nn))-min(l1_norm_c(:,nn)))^2;
    crit3 = abs(tv_penalty(:,nn)-min(tv_penalty(:,nn))).^2 /(max(tv_penalty(:,nn))-min(tv_penalty(:,nn)))^2;
    crit = crit1+crit2+(1+3*noise_std(nn))*crit3;
    select_ind(nn) = find( (crit == min(crit)),1 );
    gamma_select(nn) = P.lambda_values(select_ind(nn));
    % Plot L-curves
    subplot(3,4,nn)
    plot(crit1+crit2,...
         crit3)
    hold on
    plot(crit1(select_ind(nn))+crit2(select_ind(nn)),...
         crit3(select_ind(nn)),'or')
    plot(crit1(s_i(nn))+crit2(s_i(nn)),...
         crit3(s_i(nn)),'sb')
    [awmv_rmse(select_ind(nn),nn),...
     awmv_rmse(s_i(nn),nn)]
    
    hold off
    xlabel('MSE+l_1-norm')
    ylabel('TV Penalty')

end

   %% Show L-curves and selected parameters
nn = 3;

[~,s_i] = min(awmv_rmse);
select_ind = zeros(NN,1);
gamma_select = zeros(NN,T);
figure(111)
for nn = 1:NN
    crit1 = (0.5*mse_c(:,nn));
    crit2 = (l1_norm_c(:,nn));
    crit3 = log(tv_penalty(:,nn));
    crit = crit1+crit2+crit3;
    select_ind(nn) = find( (crit == min(crit)),1 );
    gamma_select(nn) = P.lambda_values(select_ind(nn));
    % Plot L-curves
    subplot(3,4,nn)
    plot(log(crit1+crit2),...
         crit3)
    hold on
    plot(log(crit1(select_ind(nn))+crit2(select_ind(nn))),...
         crit3(select_ind(nn)),'or')
    plot(log(crit1(s_i(nn))+crit2(s_i(nn))),...
         crit3(s_i(nn)),'sb')
    [awmv_rmse(select_ind(nn),nn),...
     awmv_rmse(s_i(nn),nn)]
    
    hold off
    xlabel('MSE+l_1-norm')
    ylabel('TV Penalty')

end
figure(444)
hold on
plot(crit1)
plot(crit2)
%% Plot the pairs of AWMV fits 
figure(222)
for nn = 1:NN
    subplot(3,4,nn)
    hold on
    plot(theta_stds1,'k')
    plot(awmv_indep(:,nn),'g')
    plot(awmv_all(select_ind(nn),:,nn),'r')
    plot(awmv_all(s_i(nn),:,nn),'b')
end
   
   
%% Plot AWMV curves for different parameter values
[param_awmv_fig, ~] = tight_subplot(5,3,[.005 .005],[.01 .01],[.01 .01]);
nn = 1
for i = 1:15
    axes(param_awmv_fig(i))
    hold on
    plot(theta_stds1,'k','Linewidth',2)
    plot(awmv_indep(:,nn),'b','Linewidth',2)
    plot(awmv_all(i,:,nn),'r','Linewidth',2)
    legend(num2str(lambda2_values(nn,i)),'Location','Northwest')
end

%%
% Plot RMSE in AWMV
figure(33)
hold on
plot(awmv_rmse_indep,'b','Linewidth',2)
plot(awmv_rmse_coupled,'r','Linewidth',2)
legend('Sparse Model','VDF Smoothed Model','Location','Best')
ylabel('Relative MSE in AWMV')
xlabel('\sigma_{noise}')

% Waterfall Data + Fits
fits = zeros(N,T,NN);
for nn = 1:NN
    dset_name = ['singlePeak_noise',num2str(nn)];
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    load(fullfile(output_dir,[dset_name,'_',num2str(s_i(nn)),'_','time',num2str(5)]))
    for time = 1:T
        x = X_hat(:,:,time);
        fits(:,time,nn) = Ax_ft_1D(A0ft_stack,x);    
    end
end

%% Waterfall Plots
figure(3)

nn = 3;
dset_name = ['singlePeak_noise',num2str(nn)];
output_name = '_coupled_CGTV';
output_subdir = [dset_name,output_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
load(fullfile(output_dir,[dset_name,'_',num2str(s_i(nn)),'_','time',num2str(5)]))

% axes(ha_waterfall(1))
subplot(3,2,1)
waterfall(B')
zlim([0 1.5])
ylabel('time')
xlabel('\eta')
zlabel('Intensity')
% axes(ha_waterfall(2))
subplot(3,2,2)
waterfall(fits(:,:,nn)')
zlim([0 1.5])
ylabel('time')
xlabel('\eta')
zlabel('Intensity')

nn = 8;
dset_name = ['singlePeak_noise',num2str(nn)];
output_name = '_coupled_CGTV';
output_subdir = [dset_name,output_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
load(fullfile(output_dir,[dset_name,'_',num2str(s_i(nn)),'_','time',num2str(5)]))

% axes(ha_waterfall(1))
subplot(3,2,3)
waterfall(B')
zlim([0 1.5])
ylabel('time')
xlabel('\eta')
zlabel('Intensity')
% axes(ha_waterfall(2))
subplot(3,2,4)
waterfall(fits(:,:,nn)')
zlim([0 1.5])
ylabel('time')
xlabel('\eta')
zlabel('Intensity')

nn = 11;
dset_name = ['singlePeak_noise',num2str(nn)];
output_name = '_coupled_CGTV';
output_subdir = [dset_name,output_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
load(fullfile(output_dir,[dset_name,'_',num2str(s_i(nn)),'_','time',num2str(5)]))

% axes(ha_waterfall(1))
subplot(3,2,5)
waterfall(B')
zlim([0 1.5])
ylabel('time')
xlabel('\eta')
zlabel('Intensity')
% axes(ha_waterfall(2))
subplot(3,2,6)
waterfall(fits(:,:,nn)')
zlim([0 1.5])
ylabel('time')
xlabel('\eta')
zlabel('Intensity')
%{
    select_ind = zeros(T,1);
    mse_indep = zeros(T,1);
    mse_coupled = zeros(T,1);
    l1_norm = zeros(M,T);
    awmv_az_indep = zeros(T,1);
    awmv_az_coupled = zeros(T,1);
    for time = 1:T
        x = X_init(:,:,time);
        fit = Ax_ft_1D(A0ft_stack,x);
        mse_indep(time) = norm(fit-B(:,time))^2;
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_az_indep(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;

        x = X_hat(:,:,time);
        fit = Ax_ft_1D(A0ft_stack,x);
        mse_coupled(time) = norm(fit-B(:,time))^2;
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_az_coupled(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end

    mse_indep_time = sqrt(sum(mse_indep));
    mse_coupled_time = sqrt(sum(mse_coupled));
    awmv_mse_indep(nn) = norm(awmv_az_indep-theta_stds1)/norm(theta_stds1);
    awmv_mse_coupled(nn) = norm(awmv_az_coupled-theta_stds1)/norm(theta_stds1);
    
    axes(ha12(nn))
    hold on
    plot(awmv_az_indep)
    plot(awmv_az_coupled)
    plot(theta_stds1)
    legend(num2str(awmv_mse_indep(nn)),num2str(awmv_mse_coupled(nn)),'\sigma',...
           'Location','Best')
end

%% Plot selected awmv

awmvs_fig = figure(11);
[ha12, ~] = tight_subplot(4,3,[.005 .005],[.01 .01],[.01 .01]);
j_select = [1,4;1,4;1,4;2,4;2,4;2,4;...
            1,2;1,2;1,2;1,2;1,2];
awmv_mse_indep = zeros(11,1);
awmv_mse_coupled = zeros(11,1);
for nn = 1
%     j = j_select(ii,1);
    jj = j_select(nn,2);
    % Load coupled solution
    dset_name = ['singlePeak_noise',num2str(nn)];
    output_name = '_coupled_CGTV';
    output_subdir = [dset_name,output_name];
    dataset =  fullfile(top_dir,dset_name);
    output_dir  = fullfile(top_dir,output_subdir);
    load(fullfile(output_dir,[dset_name,'_',num2str(num_ims),'_','all2']))
    load(fullfile(output_dir,[dset_name,'_',num2str(j),'_','time',num2str(jj)]))
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
        crit = abs(l1_norm(:,time)*0.45).^2 + abs(mse_indep(:,time)).^2;
        select_ind(time) = find( (crit == min(crit)),1 );
        lambda1_select(time) = P.lambda_values(select_ind(time));
        x = X_indep(:,:,select_ind(time),time);
        X_init(:,:,time) = x;
    end

    select_ind = zeros(T,1);
    mse_indep = zeros(T,1);
    mse_coupled = zeros(T,1);
    l1_norm = zeros(M,T);
    awmv_az_indep = zeros(T,1);
    awmv_az_coupled = zeros(T,1);
    for time = 1:T
        x = X_init(:,:,time);
        fit = Ax_ft_1D(A0ft_stack,x);
        mse_indep(time) = norm(fit-B(:,time))^2;
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_az_indep(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;

        x = X_hat(:,:,time);
        fit = Ax_ft_1D(A0ft_stack,x);
        mse_coupled(time) = norm(fit-B(:,time))^2;
        az_signal = squeeze(sum(x,1));
        var_sum = squeeze(sum(az_signal(:)));
        awmv_az_coupled(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    end

    mse_indep_time = sqrt(sum(mse_indep));
    mse_coupled_time = sqrt(sum(mse_coupled));
    awmv_mse_indep(nn) = norm(awmv_az_indep-theta_stds1)/norm(theta_stds1);
    awmv_mse_coupled(nn) = norm(awmv_az_coupled-theta_stds1)/norm(theta_stds1);
    
    axes(ha12(nn))
    hold on
    plot(awmv_az_indep)
    plot(awmv_az_coupled)
    plot(theta_stds1)
    legend(num2str(awmv_mse_indep(nn)),num2str(awmv_mse_coupled(nn)),'\sigma',...
           'Location','Best')
end

figure(2)
hold on
plot(awmv_mse_indep)
plot(awmv_mse_coupled)
legend('indep','coupled','Location','Best')
ylabel('Relative MSE in AWMV')
xlabel('\sigma_{noise}')

%% Plot L-curve
% figure(3)
% plot(l1_norm,mse_indep,'o-')
% % waterfall(1:T,l1_norm,rse_indep)
% xlabel('l_1 norm')
% ylabel('Rel Error')


%% Plot awmv, rel err 
% awmv_fig = figure(2);
% hold on
% plot(theta_stds1,theta_stds1,'-','Linewidth',2)
% % plot(theta_stds1,awmv_az_indep,'Linewidth',1)
% plot(theta_stds1,awmv_az_indep_CG,'Linewidth',1)
% plot(theta_stds1,awmv_az,'Linewidth',1)
% xlabel('\sigma')
% ylabel('AWMV')
% 
% mse1 = norm(theta_stds1-awmv_az_indep);
% mse2 = norm(theta_stds1-awmv_az_indep_CG);
% mse3 = norm(theta_stds1-awmv_az);
% 
% % mse1 = norm(theta_stds1-awmv_az_indep)/norm(theta_stds1);
% % mse2 = norm(theta_stds1-awmv_az_indep_CG)/norm(theta_stds1);
% % mse3 = norm(theta_stds1-awmv_az)/norm(theta_stds1);
% 
% legend('truth',['indep- ',num2str(mse2)],...
%                 ['coupled- ',num2str(mse3)],...
%                 'location','best')
%             
% % legend('truth',['indep- ',num2str(mse1)],...
% %                 ['indep_{CG}- ',num2str(mse2)],...
% %                 ['coupled- ',num2str(mse3)],...
% %                 'location','best')
% 
% rel_err_fig = figure(3);
% hold on
% plot(theta_stds1,rse_indep,'-','Linewidth',2)
% plot(theta_stds1,rse_indep_CG,'-','Linewidth',2)
% plot(theta_stds1,rse_coup,'-','Linewidth',1)
% legend('indep','indep_{CG}','coupled','location','best')
% %end

%% Check awmv
% close all
% ii = 5;
% i = lambda_inds(ii);
% dset_name = ['singlePeak_noise',num2str(ii)];
% % Output dirs
% output_name = '_indep_ISM';
% output_subdir = [dset_name,output_name];
% output_dir  = fullfile(top_dir,output_subdir);
% 
% load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','all']))
% awmv_fig = figure(2);
% hold on
% plot(theta_stds1,theta_stds1,'-','Linewidth',2)
% plot(theta_stds1,awmv_az,'Linewidth',1)
% xlabel('\sigma')
% ylabel('AWMV')
% awmv_fig.Position = [800 700 200 200];
%     set(gca,'XTickLabel',[]);
%     set(gca,'YTickLabel',[]);
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'XColor', 'none','YColor','none')

%% Produce 3 example fits
fits_fig = figure(1);
[ha2, ~] = tight_subplot(11,3,[.005 .005],[.01 .01],[.01 .01]); 
% lambda_inds = [15,17,19,21,...
%                22,22,22,23,...
%                23,23,24];
im_ind=1;
for nn = 1:11
j = j_select(nn,1);
jj = j_select(nn,2);
dset_name = ['singlePeak_noise',num2str(nn)];
indep_name = '_indep_ISM';
output_name = '_coupled_CGTV';
output_subdir = [dset_name,output_name];
indep_subdir = [dset_name,indep_name];
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);
indep_dir  = fullfile(top_dir,indep_subdir);
load(fullfile(indep_dir,[dset_name,'_',num2str(num_ims),'_','all2']))
load(fullfile(output_dir,[dset_name,'_',num2str(j),'_','time',num2str(jj)]))

awmv_az = zeros(T,1);
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
for t = [1,25,50]
    x = X_hat(:,:,t);
    fit = Ax_ft_1D(A0ft_stack,x);
    az_signal = squeeze(sum(x,1));
    var_sum = squeeze(sum(az_signal(:)));
    awmv_az(t) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
    
    axes(ha2(im_ind))
    hold on
    plot(B(:,t),'Linewidth',2)
    plot(fit,'Linewidth',1)
    im_ind = im_ind+1;
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'XColor', 'none','YColor','none')
end
fits_fig.Position = [800 800 300 100];
end
%}