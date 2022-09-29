function fig1 = plotDeltaAWMVbyOmega(x_time,delta_awmv_all,savePNG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
lineStyles = {'-','--',':','-.'};
rings = {'\{004\}','\{021\}','\{112\}','\{020\}'};
for om = 2:5
    fig1 = figure(om);
    hold on
    for ring = 1:4
        plot(x_time,delta_awmv_all(:,ring,om-1),...
                'LineWidth',1.5,...
                'LineStyle',lineStyles{ring})
        ylim([0,1.4])   
        ylabel('\Delta AWMV in \eta (degrees)','FontSize',14)
        xlabel('time (s)','FontSize',14)
        title(['\Delta\omega = ',num2str(om),'^\circ'],'FontSize',14)
    end
    c_leg = legend(rings,'Location','SE','FontSize',14);
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','Ring indices \{hkl\}')
    grid on
    if savePNG
        saveas(fig1,['awmvs_omega',num2str(om),'ring1234.png'])
    end
end

