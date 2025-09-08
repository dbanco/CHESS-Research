function trials_figures_subplot_sim1(fig_dir,num_trials,data_name,meanSNR,error_stats_indep,error_stats_of)

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
% xlabel('SNR (dB)','Fontsize',14)
ylabel('${\bf f}$ Error','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on
% fig1.Position(3:4) = [400 200];

% saveas(gcf, fullfile(fig_dir,[data_name,'_true_recon_error.png']));

% Dictionary Error Figure
subplot(7,1,2)
hold on
errorbar(meanSNR,error_stats_indep.avg_Derror,...
                 error_stats_indep.std_Derror,'s--')
errorbar(meanSNR,error_stats_of.avg_Derror,...
                 error_stats_of.std_Derror,'o-')
% title(title_str)
% xlabel('SNR (dB)','Fontsize',14)
ylabel('${\bf D}$ Error','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
% legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%        'MCDL-OF')
grid on
% fig2.Position(3:4) = [400 200];

% saveas(gcf, fullfile(fig_dir,[data_name,'_dict_rel_error.png']));

% Log Penalty Figure
subplot(7,1,3)
hold on
errorbar(meanSNR,error_stats_indep.avg_log_penalty,...
                 error_stats_indep.std_log_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_log_penalty,...
                 error_stats_of.std_log_penalty,'o-')
% title(title_str)
% xlabel('SNR (dB)','Fontsize',14)
ylabel('LogP','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
% legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%        'MCDL-OF')
grid on
fig1.Position(3:4) = [400 700];

% OF+HS Penalty Figure
subplot(7,1,4)
hold on
errorbar(meanSNR,error_stats_indep.avg_ofhs_penalty,...
                 error_stats_indep.std_ofhs_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_ofhs_penalty,...
                 error_stats_of.std_ofhs_penalty,'o-')
% title(title_str)
% xlabel('SNR (dB)','Fontsize',14)
ylabel('$ \Omega_{T} $','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
% legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%        'MCDL-OF')
grid on
fig1.Position(3:4) = [400 700];

% OF Penalty Figure
subplot(7,1,5)
hold on
errorbar(meanSNR,error_stats_indep.avg_of_penalty,...
                 error_stats_indep.std_of_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_of_penalty,...
                 error_stats_of.std_of_penalty,'o-')
% title(title_str)
% xlabel('SNR (dB)','Fontsize',14)
ylabel('$\Omega_{T,1}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
% legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%        'MCDL-OF')
grid on
% fig2.Position(3:4) = [400 200];

% saveas(gcf, fullfile(fig_dir,[data_name,'_of.png']));

% HS Penalty Figure
subplot(7,1,6)
hold on
errorbar(meanSNR,error_stats_indep.avg_hs_penalty,...
                 error_stats_indep.std_hs_penalty,'s--')
errorbar(meanSNR,error_stats_of.avg_hs_penalty,...
                 error_stats_of.std_hs_penalty,'o-')
% title(title_str)
ylabel('$ \Omega_{T,2}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
% legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%        'MCDL-OF')
grid on
% fig2.Position(3:4) = [400 200];

% saveas(gcf, fullfile(fig_dir,[data_name,'_hs.png']));

% X-metric Figure
subplot(7,1,7)
hold on
errorbar(meanSNR,error_stats_indep.avg_x_metric,...
                 error_stats_indep.std_x_metric,'s--')
errorbar(meanSNR,error_stats_of.avg_x_metric,...
                 error_stats_of.std_x_metric,'o-')
% title(title_str)
xlabel('SNR (dB)','Fontsize',14)
ylabel('$ {\bf X}$ NN-E','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
% legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%        'MCDL-OF')
grid on

saveas(gcf, fullfile(fig_dir,[data_name,'_subplot.png']));

end