function figOut = plotDeltaAWMVbyRing(x_time,delta_awmv_all,savePNG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
lineStyles = {'-','--',':','-.'};
rings = {'\{004\}','\{021\}','\{112\}','\{020\}'};
oms = {'2^\circ','3^\circ','4^\circ','5^\circ'};
for ring = 1:4
    figOut = figure(ring);
    hold on
    for om = 2:5
        plot(x_time,delta_awmv_all(:,ring,om-1),...
                'LineWidth',1.5,...
                'LineStyle',lineStyles{om-1})
        ylim([0,1.4])   
        ylabel('\Delta AWMV in \eta (degrees)','FontSize',14)
        xlabel('time (s)','FontSize',14)
        title(rings{ring},'FontSize',14)
    end
    c_leg = legend(oms,'Location','SE','FontSize',14);
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','\Delta\omega')
    grid on
    if savePNG
        saveas(figOut,['awmvs_ring',num2str(ring),'omega1234.png'])
    end
end

