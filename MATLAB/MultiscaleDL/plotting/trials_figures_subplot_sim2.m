function trials_figures_subplot_sim2(fig_dir,num_trials,data_name,meanSNR,error_stats_indep,error_stats_of)

title_str = ['Average over ', num2str(num_trials),' trials'];

% True Recon Error Figure
fig1 = figure;
subplot(7,1,1)
hold on
isequal(size(meanSNR), size(error_stats_indep.avg_true_error), size(error_stats_indep.std_true_error))
disp(size(meanSNR))
disp(size(error_stats_indep.avg_true_error))
disp(size(error_stats_indep.std_true_error))
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
subplot(7,1,2)
hold on
errorbar(meanSNR,error_stats_indep.avg_Derror,...
                 error_stats_indep.std_Derror,'s--')
errorbar(meanSNR,error_stats_of.avg_Derror,...
                 error_stats_of.std_Derror,'o-')

ylabel('${\bf D}$ Error','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on

% Log Penalty Figure
subplot(7,1,3)
hold on
errorbar(meanSNR,error_stats_indep.avg_log_penalty,...
                 error_stats_indep.std_log_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_log_penalty,...
                 error_stats_of.std_log_penalty,'o-')

ylabel('LogP','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on
fig1.Position(3:4) = [400 700];

% X Metric Figure
subplot(7,1,4)
hold on
errorbar(meanSNR,error_stats_indep.avg_x_metric,...
                 error_stats_indep.std_x_metric,'s--')
errorbar(meanSNR,error_stats_of.avg_x_metric,...
                 error_stats_of.std_x_metric,'o-')

ylabel('${\bf X} Metric$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on


% OF+HS Penalty Figure
subplot(7,1,5)
hold on
errorbar(meanSNR,error_stats_indep.avg_ofhs_penalty,...
                 error_stats_indep.std_ofhs_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_ofhs_penalty,...
                 error_stats_of.std_ofhs_penalty,'o-')

ylabel('$ \Omega_{T} $','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on

% OF Penalty Figure
subplot(7,1,6)
hold on
errorbar(meanSNR,error_stats_indep.avg_of_penalty,...
                 error_stats_indep.std_of_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_of_penalty,...
                 error_stats_of.std_of_penalty,'o-')

ylabel('$\Omega_{T,1}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
grid on

% HS Penalty Figure
subplot(7,1,7)
hold on
errorbar(meanSNR,error_stats_indep.avg_hs_penalty,...
                 error_stats_indep.std_hs_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_hs_penalty,...
                 error_stats_of.std_hs_penalty,'o-')

xlabel('SNR (dB)','Fontsize',14)
ylabel('$ \Omega_{T,2}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)

grid on

fig1.Position(3:4) = [400 700];

saveas(gcf, fullfile(fig_dir,[data_name,'_subplot.png']));

end