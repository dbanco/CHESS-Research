% MMPAD metadata
mmpad_dims = [396,265];
rings = {'\{004\}','\{021\}','\{112\}','\{020\}'};
% d_space = [1.27773,1.24845,1.23272,1.17138];
% two_theta = [6.89132,7.05316,7.14328,7.5179];

d_space = [1.17138,1.23272,1.24845,1.27773];
two_theta = [7.5179,7.14328,7.05316,6.89132];

start_pixel = [5,210,268,355];

total_time = 136.25;
fps = 4;
x_time = 0:1/fps:total_time;

pixel_side = 150e-6;
gap_width = 0.75e-3;

detec_dist = 4.6;
detect_angle = 14.6;

% Compute circumeference of rings
radius = 4.6*tan(pi/180*two_theta);
circum = 2*pi*radius;

% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;
% pixel_angle = 180/pi*atan(pixel_side/detec_dist);

eta_start = 14.6 - 130*pixel_angle
eta_end =   14.6 + 130*pixel_angle
top_dirs = {'D:\MMPAD_data_nr1','D:\MMPAD_data_nr2'};
data_tags = {'_nr1','_nr2'};
num_rings = [4,3];
lineStyles = {'-','--',':','-.'};
d_num = 1
top_dir = top_dirs{d_num};
data_tag = data_tags{d_num};
dset_name = 'ring%i_zero';
prefix = 'mmpad_img_%i.mat';

%% Compute center of mass
R = num_rings(d_num);
N = 546;
COM = zeros(R,N);
angleCOM = zeros(R,N);
for ring = 1:R
    for img = 1:N
        ring_dir = sprintf(dset_name,ring);
        img_file = sprintf(prefix,img);
        load(fullfile(top_dir,ring_dir,img_file))

        [n,m] = size(polar_image);
        radial_img = sum(polar_image,2);
        radial_coord = (1:n)';
        total_mass = sum(radial_img);
        COM(ring,img) = sum(radial_img.*radial_coord./total_mass);
    end 
end
close all
ring_pos = zeros(4,N);
rel_ring_pos = zeros(4,N);
for ring = 1:R
    ring_pos(ring,:) = COM(ring,:) + start_pixel(ring);
end
for ring = 1:R
    rel_ring_pos(ring,:) = ( ring_pos(4,1)-ring_pos(ring,:) )*pixel_side;
    approx_rel_ring_pos(ring,:) = rel_ring_pos(ring,:) + radius(4) - 5*pixel_side;
    angleCOM(ring,:) = atan(approx_rel_ring_pos(ring,:)./detec_dist);
end
plot(rel_ring_pos')
%% Plot all CoM and convert to angle units
close all
outName =  ['allCoM',data_tag,'.png'];
datadir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
fig_CoM = figure(9);
delta_angleCOM = zeros(R,N);
for ring = 1:R      
    hold on
    delta_angleCOM(ring,:) = angleCOM(ring,:) - angleCOM(ring,1);
    plot(x_time,180/pi*delta_angleCOM(ring,:),...
        'LineWidth',1.5,...
        'LineStyle',lineStyles{ring})
    ylabel('\Delta CoM in 2\theta (degrees)','FontSize',14)
    xlabel('time (s)','FontSize',14)

end
c_leg = legend(rings,'Location',' Best','FontSize',14);
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')
grid on
saveas(fig_CoM,fullfile(datadir,outName))

%% Plot strain
close all
outName =  ['allStrain',data_tag,'.png'];
datadir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
fig_CoM = figure(9);
strain = zeros(R,N);
for ring = 1:R
    hold on
    strain(ring,:) = -cot(angleCOM(ring,1)/2)*delta_angleCOM(ring,:)/2;
    plot(x_time,strain(ring,:),'LineWidth',1.5,'LineStyle',lineStyles{ring})
    ylabel('\Delta Strain','FontSize',14)
    xlabel('time (s)','FontSize',14)

end
c_leg = legend(rings,'Location','SouthWest','FontSize',14);
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')
grid 
%     ylim([-0.00005 0.002])
%     xlim([2.5 9])

yyaxis right
% Import load curve
loadCurveFile = 'C:\Users\dpqb1\CHESS-Research\PureTiTD_loading.xlsx';
sheet = 'PureTiTD_nr1_c';
xlsRange = 'I4:I127';
loadCurve = xlsread(loadCurveFile,sheet,xlsRange);
loadTime = 2.3:1.088757396:137.3

lCurve = plot(loadTime,loadCurve,'-o',...
    'Linewidth',2,...
    'Color','#aaa',...
    'MarkerIndices',1:10:500,...
    'MarkerSize',6);

ylim([min(loadCurve)-50 458])
%     ylim([min(loadCurve)-2 222])
ylabel('Stress (MPa)')
ax = gca;
ax.YAxis(2).Color = 'k';
ah1=axes('position',get(gca,'position'),'visible','off');
leg2 = legend(ah1,lCurve,'Stress','Location','Best','FontSize',14);

c_leg.String = rings

saveas(fig_CoM,fullfile(datadir,outName))

%% Coupled and CoM
fileBase = ['ring%i_zero_mirror_coupled_awmv.mat'];
close all
select_inds = [29 23 25 25];
allComVSawmv = figure(9);
hold on
for ring = 1:R
    load(fullfile(datadir,sprintf(fileBase,2)))
    % Rescale center of Mass
    minAWMV = min(awmv_az(select_inds(ring),:)*pixel_angle(ring));
    maxAWMV = max(awmv_az(select_inds(ring),:)*pixel_angle(ring));
    scaleCOM = (angleCOM(ring,:)-min(angleCOM(ring,:)))/(max(angleCOM(ring,:))-min(angleCOM(ring,:)));
    scaleCOM = (scaleCOM + minAWMV/maxAWMV)*maxAWMV;        
end


%% Plot all coupled delta awmv
close all
fileBase = ['ring%i_zero_mirror_coupled_awmv.mat'];
outName =  ['coupledDeltaAWMV',data_tag,'.png'];
fig_out = figure(55);
select_inds = [29 23 25 25];
delta_awmv_coupled = zeros(546,4);
for ring = 1:R

    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    awmv_plot = awmv_az(select_inds(ring),:)*pixel_angle(ring);
    delta_awmv_coupled(:,ring) = awmv_plot - min(awmv_plot);
    plot(x_time,delta_awmv_coupled(:,ring),...
        'LineWidth',1.5,...
        'LineStyle',lineStyles{ring})
    ylabel('\Delta AWMV in \eta (degrees)','FontSize',14)
    xlabel('time (s)','FontSize',14)

end
c_leg = legend(rings,'Location','Best','FontSize',14);
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')
grid on
% yyaxis right
% % Import load curve
% loadCurveFile = 'C:\Users\dpqb1\CHESS-Research\PureTiTD_loading.xlsx';
% sheet = 'PureTiTD_nr1_c';
% xlsRange = 'I4:I127';
% loadCurve = xlsread(loadCurveFile,sheet,xlsRange);
% loadTime = 2.3:1.088757396:137.3
% lCurve = plot(loadTime,loadCurve,'-o','Linewidth',2,'Color','#aaa',...
%     'MarkerIndices',1:10:500,...
%     'MarkerSize',6);
% ylim([min(loadCurve) 450])
% ylabel('Stress (MPa)')
% ax = gca;
% ax.YAxis(2).Color = 'k';
% ah1=axes('position',get(gca,'position'),'visible','off');
% leg2 = legend(ah1,lCurve,'Stress','Location','Best','FontSize',14);
% grid on
% c_leg.String = rings;
saveas(fig_out,fullfile(datadir,outName))

       %% Plot strain vs awmv
close all
outName =  ['strain_vs_awmv',data_tag,'.png'];
datadir = 'C:\Users\dpqb1\Desktop\AWMV_mirror_Figures';
fig_CoM = figure(9);

for ring = 1:R
    hold on
    plot(strain(ring,:),delta_awmv_coupled(:,ring),...
        'LineWidth',1.5,...
        'LineStyle',lineStyles{ring})
    xlabel('\Delta Strain','FontSize',14)
    ylabel('\Delta AWMV in \eta (degrees)','FontSize',14)
    xlim([-0.2e-3 2.6e-3])
    ylim([-0.1 1.3])

end
c_leg = legend(rings,'Location','West','Fontsize',14);
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')
grid on
saveas(fig_CoM,fullfile(datadir,outName))


%% Plot indep awmvs against coupled
close all
outName = ['indep_coupled_AWMV',data_tag,'.png'];
fig_out = figure(7);
colors = colororder;

for ring = 1:R

    load(fullfile(datadir,sprintf(fileBase,ring)))

%         subplot(1,4,ring)
    subplot(2,2,ring)
    fontSz = 14;
    hold on
    awmv_plot = awmv_az_init*pixel_angle(ring);
    delta_awmv = awmv_plot - min(awmv_plot);
    plot(x_time,delta_awmv,'LineWidth',1,'Color',[0.5 0.5 0.5])
    plot(x_time,delta_awmv_coupled(:,ring),'d',...
        'MarkerIndices',1:30:546,...
        'LineWidth',1.5,...
        'Color',colors(ring,:),...
        'LineStyle',lineStyles{1})
    ylabel('\Delta AWMV in \eta (degrees)','FontSize',fontSz)
    xlabel('time (s)','FontSize',fontSz)
    ylim([0,1.5])
    SE = [max(xlim) min(ylim)]-[diff(xlim)*0.25 -diff(ylim)*0.1];
    text(SE(1), SE(2), rings{ring},'FontSize',fontSz)
    legend('\gamma = 0','\gamma = \gamma*',...
       'FontSize',16,'EdgeColor',[1 1 1],'location','South')
end
legend('\gamma = 0','\gamma = \gamma*',...
       'FontSize',16,'EdgeColor',[1 1 1],'location','South')
%     fig_out.Position = [700 700 1500 270]
fig_out.Position = [500 500 700 480];
saveas(fig_out,fullfile(datadir,outName))


%% Plot all Indep AWMVs
close all
outName = ['indepDeltaAWMV',data_tag,'.png'];
fig_out = figure(7);
for ring = 1:R

    load(fullfile(datadir,sprintf(fileBase,ring)))
    hold on
    plot(x_time,awmv_az_init*pixel_angle(ring),...
        'LineWidth',1.5,...
        'LineStyle',lineStyles{ring})
    ylabel('AWMV in \eta (degrees)')
    xlabel('time (s)')

end
c_leg = legend(rings,'Location','Best');
ctitle = get(c_leg,'Title');
set(ctitle,'String','Ring indices \{hkl\}')
grid on
saveas(fig_out,fullfile(datadir,outName))
