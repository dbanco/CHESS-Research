close all

%% Load data
load('D:\CHESS_data\al7075_311_polar_fit3\fista_fit_1_75.mat')
[~,NN] = size(polar_image); 
pixel_angle = NN/360;
ind1 = 100;
ind2 = 480;
n = ind2 - ind1 + 1;
% b = 1e-4*(sum(polar_image(:,ind1:ind2),1))';
b =   1e-1*log(sum(polar_image(:,ind1:ind2),1))';
b = b - min(b(:))';
%%

% Data/Dictionary Parameters
N = numel(b);
K = 20;

% Zero padding and mask
zPad = 0;
zMask = [];

% zMask = [1:zPTad,(n+zPad+1):(n+2*zPad)];

P.dataScale = 1e-4;
P.num_theta = N;
P.basis = 'norm2';
P.cost = 'l1';
P.num_var_t = K;
P.var_theta = linspace(0.5,50,P.num_var_t).^2;

% algorithm parameters
P.params.rho1 = 1;
P.params.lambda1 = 0.0052;
P.params.tau = 1.05;
P.params.mu = 2;
P.params.adaptRho = 1;
P.params.alpha = 1.8;
P.params.stoppingCriterion = 'OBJECTIVE_VALUE';
P.params.maxIter = 800;
P.params.tolerance = 1e-8;
P.params.isNonnegative = 1;
P.params.zeroPad = zPad;
P.params.zeroMask = zMask;
P.params.plotProgress = 0;
P.params.verbose = 1;

% Construct dictionary
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);

x_init = zeros(size(A0ft_stack));
x_hat = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,P.params);

A0_stack = unshifted_basis_vector_stack_zpad(P);

%% Show data and decomposed fit
close all
ii = 1;
colors = jet(20);
[n,m] = size(x_hat);
indices = find(x_hat>1e-4);
[row,mm] = ind2sub([n,m],indices);
fit = Ax_ft_1D(A0ft_stack,x_hat);

decomp_fig = figure(1);
xlim([0 n])
hold on

dict_fig = figure(2);
xlim([0 n/3])
hold on

Lwidth = 2;

maxAll = max(fit)*1.1;

% ylim([0 maxAll])
% ii = ii + 1;

max_old = 0;
for k = 20:-1:1
    x_j = zeros(size(x_hat));
    x_j(:,k)= x_hat(:,k);
%         x_j(:,k+1)= x(:,k+1);
    fit_j = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_j),P.params.zeroMask);
    dictPeak = shift1D(A0_stack(:,k),floor(n/2));
    dictMax = max(dictPeak);
    dictPeak = dictPeak/dictMax;
    figure(1)
    plot(fit_j + max_old,'LineWidth',Lwidth,'Color',colors(k,:))
    if(k==2)
        stem_big = find( x_j(:,k) > mean(x_j(:,k))/100 );
        stem(stem_big,x_j(stem_big,k)*dictMax + max_old,...
            'BaseValue',max_old,...
            'Marker','o',...
            'MarkerSize',3,...
            'LineWidth',2,...
            'MarkerFaceColor','#ff2133',...
            'Color','#ff4444')
%         plot([1,n],max_old*[1,1],'Color','#ff4444')
    end
    
    figure(2)
    plot(dictPeak(floor(n*1/3):floor(n*2/3))*max(max(fit_j(:))*0.75,0.05) + max_old,...
        'Linewidth',2,'Color',colors(k,:))

    space = max(max(fit_j(:))*1.1,0.05);
    max_old = max_old + space;
end

figure(1)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'XColor', 'none','YColor','none')
plot(fit + max_old,'Linewidth',Lwidth,'Color','#888888')
space = max(max(fit(:))*1.1,0.05);

figure(2)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'XColor', 'none','YColor','none')
plot(fit(floor(n*1/3):floor(n*2/3)) + max_old,'Linewidth',Lwidth,'Color','#FFFFFF')
space = max(max(fit(:))*1.1,0.05);

max_old = max_old + space;

figure(1)
plot(b + max_old,'Linewidth',Lwidth,'Color','#333333')
decomp_fig.Position = [0 0 600 700];
saveas(decomp_fig,'C:\Users\dpqb1\Desktop\paper_figures\decomp_fig3.png')

figure(2)
plot(b(floor(n*1/3):floor(n*2./3)) + max_old,'Linewidth',Lwidth,'Color','#FFFFFF')
dict_fig.Position = [550 0 120 700];
saveas(dict_fig,'C:\Users\dpqb1\Desktop\paper_figures\dict_peaks3.png')

std_i = sqrt(P.var_theta);
colors = jet(P.num_var_t);
uvdf = sum(x_hat,1);
figure(3)
hold on
for i = 1:K
    bar_graph = bar(std_i(i)/pixel_angle,uvdf(i),0.35);
    set(bar_graph, 'FaceColor', colors(i,:))
end
xlabel('\sigma_i  (\circ)','FontSize',18)
ylabel('UVDF','FontSize',18)

%% Show dictionary peaks
% [n,m] = size(A0_stack);
% 
% dict_fig = figure(3);
% hold on
% maxAll = max(fit)*1.1;
% xlim([0 n/4])
% % ylim([0 maxAll])
% % ii = ii + 1;
% 
% max_y = 0;
% max_old = 0;
% for k = 20:-1:1
%     x_j = zeros(size(x_hat));
%     x_j(:,k)= x_hat(:,k);
% %         x_j(:,k+1)= x(:,k+1);
%     fit_j = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_j),P.params.zeroMask);
%     dictPeak = shift1D(A0_stack(:,k),floor(n/2));
%     space = max(max(fit_j(:))*1.1,0.05);
%     plot(dictPeak(floor(n*1.5/4):floor(n*2.5/4))*max(fit_j(:)) + max_old,'Linewidth',1,'color',colors(k,:))
%     max_old = max_old + space;
% end



%     plot(zeros(n,1),'Color','White','LineWidth',1.5)
%     set(gca,'Visible','off')

% figure(33)
% plot(b,'LineWidth',1.5,'Color',[0 0 0])
% hold on
% plot(fit,'LineWidth',1,'Color','red','LineStyle','--')
% ylim([0 max(b)])
% %     legend('Data','Fit')
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
% ii = ii + 1;
% 
