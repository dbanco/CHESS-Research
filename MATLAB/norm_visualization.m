% Norm visualization script

names = {'$l^2$ norm','$l^1$ norm','$l^{0.5}$ norm','$l^0$ norm'}
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




    