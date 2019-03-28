az_mean2 = 52:4:112;  %16
az_std_both = 1:1:15; %15

sig_i = 2:15;
dist_i = [2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5];

sig_i2 = 1:15;
dist_i2 =[1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9];

for i = dist_i2
    y = az_mean2(dist_i2)-50;
end


figure(1)
plot(sig_i2,y,'o-')
hold on 

% Linear regression
A = [sig_i2',ones(15,1)];
% A = sig_i2';
m = linsolve(A,y');
fit = A*m;
plot(sig_i2,fit,'-')
ylabel('Separation')
xlabel('\sigma')
legend('Data',sprintf('Fit: m = %1.2f, b = %1.2f',m(1),m(2)),'location','best')