% MMPAD metadata
mmpad_dims = [396,265];
rings = {'020','112','021','004'};
d_space = [1.27773,1.24845,1.23272,1.17138];
two_theta = [6.89132,7.05316,7.14328,7.5179];

pixel_side = 150e-6;
gap_width = 0.75e-3;

detec_dist = 4.6;
detect_angle = 14.6;

% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;

% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;

R = 4;

%% Load AWMV
datadir = 'E:\MMPAD_omega';
fileBase = 'ring%iomega%i_mirror_indep_awmv.mat';
outName =  'ring%iomega%i_indep_AWMV.png';
origBase = 'ring%i_zero_mirror_indep_awmv.mat';
for ring = 1:R
    fig_out = figure(ring);
    hold on
    for om = 2:4
        load(fullfile(datadir,sprintf(fileBase,ring,om)))
        % Rescale center of Mass
        minAWMV = min(awmv);
        maxAWMV = max(awmv);
        plot(awmv,'LineWidth',1.5)
    end
    load(fullfile(datadir,sprintf(origBase,ring)))
    plot(awmv,'LineWidth',1.5)
    title(sprintf('Ring %i',ring))
    legend('2^\circ','3^\circ','4^\circ','5^\circ')
%     saveas(fig_out,fullfile(datadir,sprintf(outName,ring,om)))
end

%{
%% PLot all awmv
datadir = 'C:\Users\dpqb1\Desktop\mmpad_rd_figures';
fileBase = 'ring%i_zero_couple_fit_data.mat';
outName =  'allAWMV.png';
fig_out = figure(5);
for ring = 1:R
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(awmv_az(select_ind,:),'LineWidth',1.5)
    ylabel('AWMV')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4','Location','Best')
saveas(fig_out,fullfile(datadir,outName))


%% Plot all Indep AWMVs
outName =  'indepAWMV.png';
fig_out = figure(7);
for ring = 1:R
    
    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(awmv_az_init,'LineWidth',1.5)
    ylabel('AWMV')
    xlabel('time')

end
legend('ring 1','ring 2','ring 3','ring 4','Location','Best')
saveas(fig_out,fullfile(datadir,outName))
%}