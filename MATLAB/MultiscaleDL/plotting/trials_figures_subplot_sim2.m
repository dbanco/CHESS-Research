function trials_figures_subplot_sim2(fig_dir,num_trials,data_name,meanSNR,error_stats_indep,error_stats_of)

title_str = ['Average over ', num2str(num_trials),' trials'];
fig1 = figure;

% True Recon Error Figure
subplot(5,2,1)
hold on
errorbar(meanSNR,error_stats_indep.avg_true_error,error_stats_indep.std_true_error,'s--')
errorbar(meanSNR,error_stats_of.avg_true_error,...
                 error_stats_of.std_true_error,'o-')
title(title_str)

ylabel('${\bf f}$ Error','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on

% Dictionary Error Figure
subplot(5,2,2)
hold on
errorbar(meanSNR,error_stats_indep.avg_Derror,...
                 error_stats_indep.std_Derror,'s--')
errorbar(meanSNR,error_stats_of.avg_Derror,...
                 error_stats_of.std_Derror,'o-')

ylabel('${\bf D}$ Error','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on

% Log Penalty Figure
subplot(5,2,3)
hold on
errorbar(meanSNR,error_stats_indep.avg_log_penalty,...
                 error_stats_indep.std_log_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_log_penalty,...
                 error_stats_of.std_log_penalty,'o-')

ylabel('LogP','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

% L0-norm Figure
subplot(5,2,4)
hold on
errorbar(meanSNR,error_stats_indep.avg_l0_norm,...
                 error_stats_indep.std_l0_norm,'s--')
errorbar(meanSNR,error_stats_of.avg_l0_norm,...
                 error_stats_of.std_l0_norm,'o-')

ylabel('$l_0$-norm','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

% OF+HS Penalty Figure
subplot(5,2,5)
hold on
errorbar(meanSNR,error_stats_indep.avg_ofhs_penalty,...
                 error_stats_indep.std_ofhs_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_ofhs_penalty,...
                 error_stats_of.std_ofhs_penalty,'o-')

ylabel('$ \Omega_{T} $','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on

% OF Penalty Figure
subplot(5,2,6)
hold on
errorbar(meanSNR,error_stats_indep.avg_of_penalty,...
                 error_stats_indep.std_of_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_of_penalty,...
                 error_stats_of.std_of_penalty,'o-')

ylabel('$\Omega_{T,1}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

% HS Penalty Figure
subplot(5,2,7)
hold on
errorbar(meanSNR,error_stats_indep.avg_hs_penalty,...
                 error_stats_indep.std_hs_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_hs_penalty,...
                 error_stats_of.std_hs_penalty,'o-')

ylabel('$ \Omega_{T,2}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on

% X Metric Figure
subplot(5,2,8)
hold on
errorbar(meanSNR,error_stats_indep.avg_x_metric,...
                 error_stats_indep.std_x_metric,'s--')
errorbar(meanSNR,error_stats_of.avg_x_metric,...
                 error_stats_of.std_x_metric,'o-')
xlabel('SNR (dB)','Fontsize',14)
ylabel('${\bf X}$ NN-E','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

% X Metric2 Figure
subplot(5,2,9)
hold on
errorbar(meanSNR,error_stats_indep.avg_x_metric2,...
                 error_stats_indep.std_x_metric2,'s--')
errorbar(meanSNR,error_stats_of.avg_x_metric2,...
                 error_stats_of.std_x_metric2,'o-')
xlabel('SNR (dB)','Fontsize',14)
ylabel('${\bf X}$ MSE','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

% Wass dist Figure
subplot(5,2,10)
hold on
errorbar(meanSNR,error_stats_indep.avg_wass_dist,...
                 error_stats_indep.std_wass_dist,'s--')
errorbar(meanSNR,error_stats_of.avg_wass_dist,...
                 error_stats_of.std_wass_dist,'o-')
xlabel('SNR (dB)','Fontsize',14)
ylabel('${\bf X}$ WD','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

fig1.Position(3:4) = [800 700];
saveas(gcf, fullfile(fig_dir,[data_name,'_subplot.png']));

end