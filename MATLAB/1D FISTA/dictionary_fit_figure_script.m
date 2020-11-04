%% Parameter selection
clear all
close all

% Parameter selection
disp('Setup parms')
P.set = 1;
% Parent directory
top_dir = 'E:\MMPAD_data';

% Input dirs
dset_name = 'ring1_zero';

% Indep dirs
indep_name = '_indep_ISM5';
indep_subdir = [dset_name,indep_name];
indep_dir = fullfile(top_dir,indep_subdir);

% Output dirs
output_name = '_coupled_CG_TVphi1';
output_subdir = [dset_name,output_name];

% Setup directories
dataset =  fullfile(top_dir,dset_name);
output_dir  = fullfile(top_dir,output_subdir);

% File Parameters
baseFileName = 'coupled_fit_%i.mat';

% Universal Parameters
% Ring sampling parameters
load(fullfile(output_dir,sprintf(baseFileName,1)))
lambda2_vals = P.lambda2_values;
% lambda2_vals = [logspace(-6,-4.1,15) logspace(-4,1,30)]
M = numel(lambda2_vals);


% Construct dictionary
baseFileName = 'indep_fit_%i_%i.mat';
load(fullfile(indep_dir,sprintf(baseFileName,1,1)))
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
A0_stack = unshifted_basis_vector_stack_zpad(P);

%% Show dictionary peaks
[n,m] = size(A0_stack);

dict_fig = figure(1);
[ha1, pos1] = tight_subplot(3,7,[0 0],[.02 .08],[.02 .02]); 
colors = jet(P.num_var_t);
hold on
for i = 1:P.num_var_t
    axes(ha1(i))
    dictPeak = A0_stack(:,i);
    plot(dictPeak+0.2,'Linewidth',1.5,'color',colors(i,:))
    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([0,1.4])
    xlim([-20,261])
    set(gca,'Visible','off')
end

%% Create color legend
leg_str = cell(m,1);
leg_fig = figure(11);
hold on
for i = 1:P.num_var_t
    dictPeak = A0_stack(:,i);
    plot(dictPeak+0.2,'Linewidth',1.5,'color',colors(i,:))
    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylim([0,1.4])
    xlim([-20,261])
    set(gca,'Visible','off')
    leg_str{i} = ['\sigma_{',num2str(i),'}'];
end
lgnd = legend(leg_str,'FontSize',18)
lgnd.NumColumns = 2;

%% Show full circulant dictionary
block_fig = figure(2);
[ha2, pos2] = tight_subplot(1,21,[0 0],[.02 .08],[.02 .02]); 
for i = 1:m
    dictBlock = zeros(n,n);
    for j = 1:n
        
        dictBlock(:,j)= shift1D(A0_stack(:,i),j);
    end
    axes(ha2(i))
    set(ha2(i), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    imagesc(dictBlock)
    colormap(ha2(i),flipud(gray.*(1-colors(i,:))+colors(i,:)) )
    set(gca,'Visible','off')
end

%% Show data and decomposed fit
load(fullfile(indep_dir,'select_indices.mat'))
dict_fig = figure(3);
[ha3, pos3] = tight_subplot(5,3,[0 0],[.02 .08],[.02 .02]); 
ii = 1;
for im_num = [1 30 40 60 100]
    param_num = select_indices(im_num);
    x_data = load( fullfile(indep_dir,sprintf(baseFileName,param_num,im_num)) );
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(x_data.P);
    [n,m] = size(A0ft_stack);

    % Fit objective
    x = x_data.x_hat;
    indices = find(x>1e-3);
    [row,mm] = ind2sub([n,m],indices);
    % fit = shift1D(Ax_ft_1D(A0ft_stack,x),131);
    fit = Ax_ft_1D(A0ft_stack,x);

    load(fullfile(dataset,[P.prefix,'_',num2str(im_num),'.mat']) )
    polar_vector = sum(polar_image,1);
    b = P.dataScale*polar_vector';
    
    axes(ha3(ii))
    max_y = 0;
    for k = 1:numel(indices)
        x_j = zeros(size(x));
        x_val = x(indices(k));
        x_j(indices(k)) = x(indices(k));
        fit_j = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_j),P.params.zeroMask);

        figure(3)
        plot(fit_j,'LineWidth',1,'Color',[colors(mm(k),:),0.5])
        hold on
        xlim([0 250])
        if max(fit_j) > max_y
            max_y = max(fit_j);
        end
        ylim([0 max_y*1.1])
    end
    plot(zeros(n,1),'Color','White','LineWidth',1.5)
    set(gca,'Visible','off')
    ii = ii + 1;

    axes(ha3(ii))
    plot(b,'LineWidth',1.5,'Color',[0 0 0])
    hold on
    plot(fit,'LineWidth',1,'Color','red','LineStyle','--')
    ylim([0 max(b)])
    legend('Data','Fit')
    set(gca,'Visible','off')
    ii = ii + 1;
end

%% Vdf figure
vdf_fig = figure(4);
[ha4, pos4] = tight_subplot(5,1,[0.05 0.05],[.05 0.05],[.1 .05]); 
ii = 1;
for im_num = [1 30 40 60 100]
    param_num = select_indices(im_num);
    x_data = load( fullfile(indep_dir,sprintf(baseFileName,param_num,im_num)) );
    A0ft_stack = unshifted_basis_vector_ft_stack_zpad(x_data.P);
    [n,m] = size(A0ft_stack);

    % Fit objective
    x = x_data.x_hat;
    indices = find(x>1e-3);
    [row,mm] = ind2sub([n,m],indices);
    % fit = shift1D(Ax_ft_1D(A0ft_stack,x),131);
    fit = Ax_ft_1D(A0ft_stack,x);

    load(fullfile(dataset,[P.prefix,'_',num2str(im_num),'.mat']) )
    polar_vector = sum(polar_image,1);
    b = P.dataScale*polar_vector';

    labs = 1:21;
    summed_x = sum(x,1);
    axes(ha4(ii))
    for bb = 1:m
        bar_graph = bar(bb,summed_x(bb));
        hold on
        set(bar_graph, 'FaceColor', colors(bb,:))
    end
    set(gca,'box','off')
    xticks(1:2:21)
%     set(gca,'Visible','off')
    ii = ii + 1;
end