% Test poisson data SNR
figure(1)
alphas = [1000 100 10 0.1 0.01 0.001];
sim = 'linear';
N = 101;
for a = 1:numel(alphas)
    [Bn,B,theta_stds,rel_err,rms] = genSimDataPoisson(N,1,alphas(a),sim);
    B = B*500*rms*3/alphas(a);
    Bn = Bn*500*rms*3/alphas(a);
 
    subplot(numel(alphas),1,a)
    hold on
    rmsS = sqrt(mean(B.^2));
    rmsN = sqrt(mean((Bn-B).^2));
    
    title([num2str( mean(sqrt(3*B)) ),'-',num2str( rmsS/rmsN )])
    plot(Bn)

end