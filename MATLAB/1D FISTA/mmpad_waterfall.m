%% Load data
clear all 
close all

rings = {'\{004\}','\{021\}','\{112\}','\{020\}'};
d_space = [1.17138,1.23272,1.24845,1.27773];
two_theta = [7.5179,7.14328,7.05316,6.89132];

% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;
pixel_side = 150e-6;
gap_width = 0.75e-3;
detec_dist = 4.6;
detect_angle = 14.6;

% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;

prefix = 'mmpad_img';

% Parent directory
ring_num = 4;
top_dir = ['D:\MMPAD_data_nr1\ring',num2str(ring_num),'_zero'];
fig_dir = 'C:\Users\dpqb1\Desktop\paper_figures\';
max_val = [2e4 2e4];

T = 546;
load(fullfile(top_dir,[prefix,'_',num2str(1),'.mat']) )
[m,n] = size(polar_image);

B = zeros(n,T);
for t = 1:T
    load(fullfile(top_dir,[prefix,'_',num2str(t),'.mat']) )
    b = squeeze(sum(polar_image,1));
    B(:,t) = b;
end
% 
% figure(1)
% waterfall((B(:,1:40)'))
% figure(2)
% waterfall((B(:,41:100)'))
% figure(3)
% waterfall((B(:,101:200)'))
% figure(4)
% waterfall((B(:,201:546)'))

%% PLot waterfall
% MMPAD metadata
mmpad_dims = [396,265];

% time axis
total_time = 136.25;
fps = 4;
x_time = 0:1/fps:total_time;

% eta axis
pixel_side = 150e-6;
gap_width = 0.75e-3;
detec_dist = 4.6;
detect_angle = 14.6;
% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;
% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;

eta_range = linspace(-261/2,261/2,261)*pixel_angle(1) + detect_angle;


clim1 = 0.5e4;
clim2 = 2e5;

figure(5)
waterfall((eta_range),(x_time(1:100)),(B(:,1:100)'))
title(rings{ring_num})
ylabel('time(s)','FontSize',16)
xlabel('\eta (\circ)','FontSize',18)
zlabel('log(Intensity)','FontSize',16)
set(gca, 'ZScale', 'log');
set(gca, 'ColorScale', 'log');
% set(gca,'CLim', [clim1 clim2])

figure(6)
waterfall((eta_range),(x_time),B')
title(rings{ring_num})
ylabel('time(s)','FontSize',16)
xlabel('\eta (\circ)','FontSize',18)
zlabel('log(Intensity)','FontSize',16)
set(gca, 'ZScale', 'log');
set(gca, 'ColorScale', 'log');
% set(gca,'CLim', [clim1 clim2])
% 
% %% Create color legend
% leg_str = cell(m,1);
% leg_fig = figure(11);
% hold on
% for i = 1:P.num_var_t
%     dictPeak = A0_stack(:,i);
%     plot(dictPeak+0.2,'Linewidth',1.5,'color',colors(i,:))
%     
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%     ylim([0,1.4])
%     xlim([-20,261])
%     set(gca,'Visible','off')
%     leg_str{i} = ['\sigma_{',num2str(i),'}'];
% end
% lgnd = legend(leg_str,'FontSize',18)
% lgnd.NumColumns = 2;
% 
% %% Show full circulant dictionary
% block_fig = figure(2);
% [ha2, pos2] = tight_subplot(1,21,[0 0],[.02 .08],[.02 .02]); 
% for i = 1:m
%     dictBlock = zeros(n,n);
%     for j = 1:n
%         
%         dictBlock(:,j)= shift1D(A0_stack(:,i),j);
%     end
%     axes(ha2(i))
%     set(ha2(i), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
%     imagesc(dictBlock)
%     colormap(ha2(i),flipud(gray.*(1-colors(i,:))+colors(i,:)) )
%     set(gca,'Visible','off')
% end
% 
% %% Show data and decomposed fit
% load(fullfile(indep_dir,'select_indices.mat'))
% dict_fig = figure(3);
% [ha3, pos3] = tight_subplot(5,3,[0 0],[.02 .08],[.02 .02]); 
% ii = 1;
% for im_num = [1 30 40 60 100]
%     param_num = select_indices(im_num);
%     x_data = load( fullfile(indep_dir,sprintf(baseFileName,param_num,im_num)) );
%     A0ft_stack = unshifted_basis_vector_ft_stack_zpad(x_data.P);
%     [n,m] = size(A0ft_stack);
% 
%     % Fit objective
%     x = x_data.x_hat;
%     indices = find(x>1e-3);
%     [row,mm] = ind2sub([n,m],indices);
%     % fit = shift1D(Ax_ft_1D(A0ft_stack,x),131);
%     fit = Ax_ft_1D(A0ft_stack,x);
% 
%     load(fullfile(dataset,[P.prefix,'_',num2str(im_num),'.mat']) )
%     polar_vector = sum(polar_image,1);
%     b = P.dataScale*polar_vector';
%     
%     axes(ha3(ii))
%     max_y = 0;
%     for k = 1:numel(indices)
%         x_j = zeros(size(x));
%         x_val = x(indices(k));
%         x_j(indices(k)) = x(indices(k));
%         fit_j = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_j),P.params.zeroMask);
% 
%         figure(3)
%         plot(fit_j,'LineWidth',1,'Color',[colors(mm(k),:),0.5])
%         hold on
%         xlim([0 250])
%         if max(fit_j) > max_y
%             max_y = max(fit_j);
%         end
%         ylim([0 max_y*1.1])
%     end
%     plot(zeros(n,1),'Color','White','LineWidth',1.5)
%     set(gca,'Visible','off')
%     ii = ii + 1;
% 
%     axes(ha3(ii))
%     plot(b,'LineWidth',1.5,'Color',[0 0 0])
%     hold on
%     plot(fit,'LineWidth',1,'Color','red','LineStyle','--')
%     ylim([0 max(b)])
%     legend('Data','Fit')
%     set(gca,'Visible','off')
%     ii = ii + 1;
% end
% 
% %% Vdf figure
% vdf_fig = figure(4);
% [ha4, pos4] = tight_subplot(5,1,[0.05 0.05],[.05 0.05],[.1 .05]); 
% ii = 1;
% for im_num = [1 30 40 60 100]
%     param_num = select_indices(im_num);
%     x_data = load( fullfile(indep_dir,sprintf(baseFileName,param_num,im_num)) );
%     A0ft_stack = unshifted_basis_vector_ft_stack_zpad(x_data.P);
%     [n,m] = size(A0ft_stack);
% 
%     % Fit objective
%     x = x_data.x_hat;
%     indices = find(x>1e-3);
%     [row,mm] = ind2sub([n,m],indices);
%     % fit = shift1D(Ax_ft_1D(A0ft_stack,x),131);
%     fit = Ax_ft_1D(A0ft_stack,x);
% 
%     load(fullfile(dataset,[P.prefix,'_',num2str(im_num),'.mat']) )
%     polar_vector = sum(polar_image,1);
%     b = P.dataScale*polar_vector';
% 
%     labs = 1:21;
%     summed_x = sum(x,1);
%     axes(ha4(ii))
%     for bb = 1:m
%         bar_graph = bar(bb,summed_x(bb));
%         hold on
%         set(bar_graph, 'FaceColor', colors(bb,:))
%     end
%     set(gca,'box','off')
%     xticks(1:2:21)
% %     set(gca,'Visible','off')
%     ii = ii + 1;
% end