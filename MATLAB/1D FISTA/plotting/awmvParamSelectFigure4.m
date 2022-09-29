function figOut = awmvParamSelectFigure4(awmv_indep,awmv_coupled,theta_stds,levels)

[~,NN] = size(awmv_indep);

figOut = figure;
figOut.Position = [1000,500,800 240];
rows = 2;
cols = 4;
[ha_awmv, ~] = tight_subplot(rows,cols,[.005 .005],[.01 .01],[.01 .01]);
for nn = 1:NN
    axes(ha_awmv(nn))
    hold on
    plot(theta_stds,'k--','Linewidth',2)
    plot(awmv_indep(:,nn),'b-','Linewidth',2)
    plot(awmv_coupled(:,nn),'r:','Linewidth',3)
    NW = [min(xlim) max(ylim)]+[diff(xlim)*0.05 -diff(ylim)*0.1];
    text(NW(1), NW(2), ['SNR=',sprintf('%1.3g',levels(nn))],...
        'LineWidth',0.1,'LineStyle','none','FontSize',14)
end
for nn = 1:rows*cols
    axes(ha_awmv(nn))
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
% Use Final plot for legend
axes(ha_awmv(NN+1))
hold on
plot(0,'k--','Linewidth',2)
plot(1,'b-','Linewidth',2)
plot(2,'r:','Linewidth',3)
legend('Truth','\gamma = 0','\gamma = \gamma*',...
       'FontSize',16,'EdgeColor',[1 1 1],'location','Northwest')

end

