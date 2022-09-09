top_dir = 'E:\MMPAD_omega\';
rings = {'ring1','ring2','ring3','ring4'};
oms = {'omega2','omega3','omega4','omega5'};

%% Plot all Indep AWMVs
 indep_dir = 'omega_mirror_indep_awmv';
for i = 1:4
    fig_out = figure(i);
    hold on
    for j = 1:4
        indepFile = [rings{i},oms{j},'_mirror_indep_awmv.mat'];
        load(fullfile(top_dir,indep_dir,indepFile))
        outName = 'indepMirrorAWMVomega.png';
        plot(awmv,'LineWidth',1.5)          
    end
    title(rings{i})
    xlabel('time (i)')
    ylabel('AWMV in \eta ')
    c_leg = legend(rings,'Location','Best');
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','Ring indices \{hkl\}')
    grid on
end

% saveas(fig_out,fullfile(datadir,outName))

%% Plot all Coupled AWMVs
coupled_dir = 'omega_mirror_coupled_awmv1';
for i = 1:4
    fig_out = figure(i);
    hold on
    for j = 1:4
        indepFile = [rings{i},oms{j},'_mirror_coupled_awmv.mat'];
        load(fullfile(top_dir,coupled_dir,indepFile))
        outName = 'coupledMirrorAWMVomega.png';
        plot(awmv_az(30,:),'LineWidth',1.5)          
    end
    title(rings{i})
    xlabel('time (i)')
    ylabel('AWMV in \eta ')
    c_leg = legend(oms,'Location','Best');
    ctitle = get(c_leg,'Title');
    set(ctitle,'String','Ring indices \{hkl\}')
    grid on
end