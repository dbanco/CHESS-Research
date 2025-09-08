function trials_figures(fig_dir,num_trials,data_name,meanSNR,error_stats_indep,error_stats_of)

title_str = ['Average over ', num2str(num_trials),' trials'];

% True Recon Error Figure
fig1 = figure;
hold on
isequal(size(meanSNR), size(error_stats_indep.avg_true_error), size(error_stats_indep.std_true_error))
disp(size(meanSNR))
disp(size(error_stats_indep.avg_true_error))
disp(size(error_stats_indep.std_true_error))
errorbar(meanSNR,error_stats_indep.avg_true_error,error_stats_indep.std_true_error,'s-')
errorbar(meanSNR,error_stats_of.avg_true_error,...
                 error_stats_of.std_true_error,'s-')
title(title_str)
xlabel('SNR (dB)','Fontsize',14)
ylabel('$\|\hat{{\bf b}}-{\bf f}\|_2$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on
fig1.Position(3:4) = [400 200];

saveas(gcf, fullfile(fig_dir,[data_name,'_true_recon_error.png']));

% Dictionary Error Figure
fig2 = figure;
hold on
errorbar(meanSNR,error_stats_indep.avg_Derror,...
                 error_stats_indep.std_Derror,'s-')
errorbar(meanSNR,error_stats_of.avg_Derror,...
                 error_stats_of.std_Derror,'s-')
title(title_str)
xlabel('SNR (dB)','Fontsize',14)
ylabel('$\|\hat{{\bf D}}-{\bf D}\|_2/\|{\bf D}\|_2$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on
fig2.Position(3:4) = [400 200];

saveas(gcf, fullfile(fig_dir,[data_name,'_dict_rel_error.png']));

% OF Penalty Figure
fig2 = figure;
hold on
errorbar(meanSNR,error_stats_indep.avg_of_penalty,...
                 error_stats_indep.std_of_penalty,'s-')
errorbar(meanSNR,error_stats_of.avg_of_penalty,...
                 error_stats_of.std_of_penalty,'s-')
title(title_str)
xlabel('SNR (dB)','Fontsize',14)
ylabel('$\Omega_{T,1}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on
fig2.Position(3:4) = [400 200];

saveas(gcf, fullfile(fig_dir,[data_name,'_of.png']));

% HS Penalty Figure
fig2 = figure;
hold on
errorbar(meanSNR,error_stats_indep.avg_hs_penalty,...
                 error_stats_indep.std_hs_penalty,'s-')
errorbar(meanSNR,error_stats_of.avg_hs_penalty,...
                 error_stats_of.std_hs_penalty,'s-')
title(title_str)
xlabel('SNR (dB)','Fontsize',14)
ylabel('$ \Omega_{T,2}$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on
fig2.Position(3:4) = [400 200];

saveas(gcf, fullfile(fig_dir,[data_name,'_hs.png']));

% OF+HS Penalty Figure
fig2 = figure;
hold on
errorbar(meanSNR,error_stats_indep.avg_ofhs_penalty,...
                 error_stats_indep.std_ofhs_penalty,'s-')
errorbar(meanSNR,error_stats_of.avg_ofhs_penalty,...
                 error_stats_of.std_ofhs_penalty,'s-')
title(title_str)
xlabel('SNR (dB)','Fontsize',14)
ylabel('$ \Omega_{T} $','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
grid on
fig2.Position(3:4) = [400 200];

saveas(gcf, fullfile(fig_dir,[data_name,'_ofhs.png']));
end