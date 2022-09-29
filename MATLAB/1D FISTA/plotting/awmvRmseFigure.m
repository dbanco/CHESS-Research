function figOut = awmvRmseFigure(snrs,awmv_rse_coup, awmv_rse_indep)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figOut = figure;
avg_awmv_coup =  squeeze(mean(awmv_rse_coup,2));
avg_awmv_indep = squeeze(mean(awmv_rse_indep,2));

std_awmv_coup =  squeeze(std(awmv_rse_coup,1,2));
std_awmv_indep = squeeze(std(awmv_rse_indep,1,2));

hold on
errorbar(snrs,avg_awmv_indep,std_awmv_indep,...
    'LineWidth',2,...
    'color','blue',...
    'LineStyle','-')
errorbar(snrs,avg_awmv_coup,std_awmv_coup,...
    'LineWidth',2,...
    'color','red',...
    'Linestyle','--')
xlabel('SNR','FontSize',16)
ylabel('Average AWMV RSE (50 trials)','FontSize',16)
legend('\gamma = 0','\gamma = \gamma*','Location','NorthEast','FontSize',16)
set(gca,'xscale','log')
grid on
figOut.Position = [1000,500,700 350];

end

