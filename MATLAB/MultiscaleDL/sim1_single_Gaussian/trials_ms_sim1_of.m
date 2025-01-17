lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];

sigmas = 0:0.01:0.1;

dataset = 'steps_matched';
penalty = 'log';
testType = 'Dflat0_Xzeros0';
topDir_of = 'E:\Outputs_sim1_of_trials\';
num_trials = 10;
fprintf('MCDL-OF Results \n')
[true_error_of,obj_array_of,data_error_of,l1_norm_of,l0_norm_of,log_penalty_of,of_penalty_of,hs_penalty_of] = sim1_trials_results(sigmas,dataset,...
                                                    penalty,testType,...
                                                    topDir_of,num_trials);
topDir_s = 'E:\Outputs_sim1_trials\';
num_trials = 10;
fprintf('MCDL Results \n')
[true_error,obj_array,data_error,l1_norm,l0_norm,log_penalty,of_penalty,hs_penalty] = sim1_trials_results(sigmas,dataset,...
                                                    penalty,testType,...
                                                    topDir_s,num_trials);
[meanSNR,noiseError] = computeSNR_noiseError(dataset);


%% True and Data Error Averaged

x = meanSNR(2:9);
y1 = noiseError(2:9);
blue = [0 0.4470 0.7410];
red =  [0.8500    0.3250    0.0980];

figure(1)
hold on
mu_s = mean(true_error(2:9,:),2);
sig_s = std(true_error(2:9,:)');
mu_of = mean(true_error_of(2:9,:),2);
sig_of = std(true_error_of(2:9,:)');

errorbar(x,mu_s,sig_s,'Color',blue)
errorbar(x,mu_of,sig_of,'Color',red)
plot(x,mu_s,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',4)
plot(x,mu_of,'s','Color',red,'MarkerFaceColor',red,'MarkerSize',4)
% plot(x,y1,'-o','MarkerSize',4)

ylabel('$\frac{1}{2}||{\bf f} - \hat{\bf f} ||_2$','interpreter','latex',...
    'FontSize',16)
xlabel('SNR','FontSize',16)
legend('MCDL','MCDL-OF')

% figure(2)
% hold on
% mu_s = mean(of_penalty(2:9,:)+hs_penalty(2:9,:),2);
% sig_s = std((of_penalty(2:9,:)+hs_penalty(2:9,:))');
% mu_of = mean(of_penalty_of(2:9,:)+hs_penalty_of(2:9,:),2);
% sig_of = std((of_penalty_of(2:9,:)+hs_penalty_of(2:9,:))');
% 
% errorbar(x,mu_s,sig_s,'Color',blue)
% errorbar(x,mu_of,sig_of,'Color',red)
% plot(x,mu_s,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',4)
% plot(x,mu_of,'s','Color',red,'MarkerFaceColor',red,'MarkerSize',4)
% % plot(x,y1,'-o','MarkerSize',4)
% 
% ylabel('OF+HS Penalty','interpreter','latex',...
%     'FontSize',16)
% xlabel('SNR','FontSize',16)
% legend('MCDL','MCDL-OF')

