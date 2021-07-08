close all

%% Load data
load('D:\CHESS_data\al7075_311_polar_fit3\fista_fit_1_75.mat')
ind1 = 190;
ind2 = 840;
n = ind2 - ind1 + 1;
% b = 1e-4*(sum(polar_image(:,ind1:ind2),1))';
b =   1e-1*log(sum(polar_image(:,ind1:ind2),1))';
b = b - min(b(:))';
figure(1)
plot(b)

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
P.params.lambda1 = 0.0002;
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

%% Show data and decomposed fit
dict_fig = figure(2);
[ha3, pos3] = tight_subplot(21,1,[0.02 0.02],[.02 .02],[.02 .02]); 
ii = 1;
colors =jet(20);
[n,m] = size(x_hat);
indices = find(x_hat>1e-4);
[row,mm] = ind2sub([n,m],indices);
fit = Ax_ft_1D(A0ft_stack,x_hat);

figure(1)
hold on
Lwidth = 2;
plot(b,'Linewidth',Lwidth)
plot(fit,'Linewidth',Lwidth-1)
maxAll = max(fit)*1.1;
xlim([0 n])
ylim([0 maxAll])
% ii = ii + 1;

max_y = 0;
for k = 1:20
    axes(ha3(ii))
    x_j = zeros(size(x_hat));
    x_j(:,k)= x_hat(:,k);
%         x_j(:,k+1)= x(:,k+1);
    fit_j = forceMaskToZero(Ax_ft_1D(A0ft_stack,x_j),P.params.zeroMask);
    
%     if max(fit_j) > 1e-3
        plot(fit_j,'LineWidth',Lwidth,'Color',[colors(k,:),0.5])

        hold on
        if max(fit_j) > max_y
            max_y = max(fit_j);
        end
    %         ylim([0 max(b)/5])
       %         ylim([0 maxAll/2]) 
        ylim([0 max(fit_j)*1.1])
        xlim([0 n])

        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        aa = gca;
        boxPosition = aa.Position;
        annotation('textbox', boxPosition, 'String', sprintf('rescaled by %2.0f',maxAll/max(fit_j)) )
        ii = ii + 1;
%     end
end
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
