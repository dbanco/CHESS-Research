%% Compare Gaussian noise data to Poisson data
[y1n, y1] = gaussian_pwlinear_2to10_problem8_gausnoise(0.04);
[y2n, y2] = gaussian_pwlinear_2to10_problem8_poisson(50);

t = 40;
figure;
hold on
plot(y1n(:,t))
plot(y2n(:,t))
plot(y1(:,t),'Linewidth',2)
% plot(y2(:,t))
legend('gaus','pois','true')

snr = norm( y1(:) )/norm( y1(:)-y1n(:) );
snr = norm( y1(:) )/norm( y1(:)-y2n(:) );