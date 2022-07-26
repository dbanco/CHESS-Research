%% Parameter selection
close all

% Parent directory
% top_dir = 'E:\PureTiRD_nr2_c_x39858';
top_dir = 'D:\Simulations';
% top_dir = '/cluster/shared/dbanco02';
mkdir(top_dir)

% Output dirs
sim_name = 'linearPoissonNoise';
alg_name = 'IndepISM';
file_name = 'lineSearchNoise';
output_dir  = fullfile(top_dir,sim_name,alg_name);
mkdir(output_dir)

a_factor = [2,5,10,20,50];

for nn = 5:numel(a_factor)
    % Setup directories
    f_name =  [file_name,'_',num2str(nn),'.mat'];
        
    % Data/Dictionary Parameters
    P.dataset = dataset;
    [P,N,K,M,T] = definePoissonP(P);
    
    % Construct dictionary
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

    % Generate data
    [B,theta_stds] = genLinearPoisson(N,T,a_factor(nn));

    %% Independent Solution
    x_init = zeros(N,K);
    X_indep = zeros(N,K,M,T);
    for i = 1:M
        P.set = i;
        P.params.lambda1 = P.lambda_values(i);
        for t = 1:T
            % Solve
            [x_hat,~,~,~,~] = convADMM_LASSO_Sherman_1D(A0ft_stack,B(:,t),x_init,P.params);  
            X_indep(:,:,i,t) = x_hat;
        end
    end
    save(fullfile(output_dir,f_name),'B','X_indep','P');
end


%% Plot fit
nn = 4;
f_name =  [file_name,'_',num2str(nn),'.mat'];
load(fullfile(output_dir,f_name))

close all
fits_fig_indep = figure(11);
 [ha11, ~] = tight_subplot(10,5,[.005 .005],[.01 .01],[.01 .01]);
awmv_az_indep = zeros(T,1);
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
%         axes(ha11(time))
%         hold on
%         plot(B(:,time))
%         plot(fit)
    end
end

%% Plot selected awmv
close all
select_ind = zeros(T,1);
for time = 1:T
    crit = abs(l1_norm(:,time)*0.5).^2 + abs(mse_indep(:,time)).^2;
    select_ind(time) = find( (crit == min(crit)),1 );
    x = X_indep(:,:,select_ind(time),time);
    az_signal = squeeze(sum(x,1));
    var_sum = squeeze(sum(az_signal(:)));
    awmv_az_indep(time) = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum;
end
select_ind
figure(1)
hold on
plot(awmv_az_indep)
plot(theta_stds)
legend('awmv','truth')

%% Plot L-curve
figure(3)
plot(l1_norm,mse_indep,'o-')
% waterfall(1:T,l1_norm,rse_indep)
xlabel('l_1 norm')
ylabel('Rel Error')

%% Plot L-curve
nn = 1
figure(4)
plot(l1_norm(15:50,nn),mse_indep(15:50,nn),'o-')
% waterfall(1:T,l1_norm,rse_indep)
xlabel('l_1 norm')
ylabel('Rel Error')


%% Plot awmv, rel err 
awmv_fig = figure(2);
hold on
plot(theta_stds,theta_stds,'-','Linewidth',2)
% plot(theta_stds1,awmv_az_indep,'Linewidth',1)
plot(theta_stds,awmv_az_indep_CG,'Linewidth',1)
plot(theta_stds,awmv_az,'Linewidth',1)
xlabel('\sigma')
ylabel('AWMV')

mse1 = norm(theta_stds1-awmv_az_indep);
mse2 = norm(theta_stds1-awmv_az_indep_CG);
mse3 = norm(theta_stds1-awmv_az);

% mse1 = norm(theta_stds1-awmv_az_indep)/norm(theta_stds1);
% mse2 = norm(theta_stds1-awmv_az_indep_CG)/norm(theta_stds1);
% mse3 = norm(theta_stds1-awmv_az)/norm(theta_stds1);

legend('truth',['indep- ',num2str(mse2)],...
                ['coupled- ',num2str(mse3)],...
                'location','best')
            
% legend('truth',['indep- ',num2str(mse1)],...
%                 ['indep_{CG}- ',num2str(mse2)],...
%                 ['coupled- ',num2str(mse3)],...
%                 'location','best')

rel_err_fig = figure(3);
hold on
plot(theta_stds1,rse_indep,'-','Linewidth',2)
plot(theta_stds1,rse_indep_CG,'-','Linewidth',2)
plot(theta_stds1,rse_coup,'-','Linewidth',1)
legend('indep','indep_{CG}','coupled','location','best')


% Check awmv
close all
ii = 5;
i = lambda_inds(ii);
dset_name = ['singlePeak_noise',num2str(ii)];
% Output dirs
output_name = '_indep_ISM';
output_subdir = [dset_name,output_name];
output_dir  = fullfile(top_dir,output_subdir);

load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','all']))
awmv_fig = figure(2);
hold on
plot(theta_stds1,theta_stds1,'-','Linewidth',2)
plot(theta_stds1,awmv_az,'Linewidth',1)
xlabel('\sigma')
ylabel('AWMV')
awmv_fig.Position = [800 700 200 200];
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'XColor', 'none','YColor','none')

% Produce 3 example fits
fits_fig = figure(1);
[ha2, ~] = tight_subplot(1,3,[.005 .005],[.01 .01],[.01 .01]); 
lambda_inds = [15,17,19,21,...
               22,22,22,23,...
               23,23,24];
ii = 3;
i = lambda_inds(ii);
dset_name = ['singlePeak_noise',num2str(ii)];
% Output dirs
output_name = '_indep_ISM';
output_subdir = [dset_name,output_name];
output_dir  = fullfile(top_dir,output_subdir);
awmv_az = zeros(T,1);
load(fullfile(output_dir,[dset_name,'_',num2str(i),'_','all']))
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
im_ind=1;
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

%}