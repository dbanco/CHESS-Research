%% Norm visualization script
names = {'$l_2$-norm','$l_1$-norm','$l_{0.5}$-norm','$l_0$-norm'};
NN = [2,1,0.5,0];
x = linspace(-1,1,501);
figure(1)
hold on
for i = 1:4
    switch i
        case 1
            y = sqrt(1-x.^2);
        case 2
            y = 1-abs(x);
        case 3
            y = (1-abs(x).^0.5).^2;
        case 4
            y = (x==0);
    end
    plot([x,x],[y,-y],'Linewidth',2)
end
legend(names,'Interpreter','latex','Fontsize',12)



%%
names = {'$l_1$-norm','$l_0$-norm'};
NN = [2,1,0.5,0];


figure(2)
hold on
plot(1,1)
plot(1,1)

c = 0.5;

d = 0.9*c;
yel_col = [0.8500 0.3250 0.0980];
yel = [1 1 1];
theta = linspace(0.4,1,4);
for i = 1:4
    y_col = yel_col*(1-theta(i)) + theta(i)*yel;
    x = linspace(-d,d,501);
    y = d-abs(x);
    plot([x,x],[y,-y],'Linewidth',2,'Color',y_col)
    d = d*0.8;
end

x = linspace(-c,c,501);
y = sqrt(c*0.4-x.^2);
plot([x,x],[y,-y],'Linewidth',2,'Color',[0 0.4470 0.7410])

x = linspace(-c,c,501);
y = c-abs(x);
plot([x,x],[y,-y],'Linewidth',2,'Color',[0.8500 0.3250 0.0980])

y = c*(x==0);
plot([x,x],[y,-y],'Linewidth',2,'Color',[0.4940 0.1840 0.5560])

x = linspace(-1,1,501);

slopes = [2 10 50];
for i = 1
    a = slopes(i);
    b = 0.5;
    y1 = x/a+b;
    plot(x,y1,'k','Linewidth',2)
end

plot(-0.2,0.4,'o','MarkerSize',25,'LineWidth',2,'Color',[0 0.4470 0.7410])
plot(0,0.5,'o','MarkerSize',25,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot(0,0.5,'o','MarkerSize',21,'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(0,0.5,'o','MarkerSize',17,'LineWidth',2,'Color',[0.4940 0.1840 0.5560])


% x = 0.2;  % X-coordinate (left)
% y = 0.8;   % Y-coordinate (top)
% width = 0.2;  % Width
% height = 0.05; % Height
% 
% annotation('textbox', [x, y, width, height], 'String', 'l_2-norm', ...
%     'FontSize', 12,'Interpreter','tex','EdgeColor', 'none');

xlim([-1.2*c,1.2*c])
ylim([-1.2*c,1.2*c])
    