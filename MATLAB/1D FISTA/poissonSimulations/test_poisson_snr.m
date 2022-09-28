% Test poisson data SNR
figure(1)
alphas = [ 1 1.1 1.25 2 3 5 10];
sim = 'linear';
N = 101;
for a = 1:numel(alphas)
    [Bn,B,theta_stds,rel_err,rms] = genSimDataPoisson2(N,1,alphas(a),sim);
    B = B*alphas(a);
    Bn = Bn*alphas(a);
 
    subplot(numel(alphas),1,a)
    hold on
    rmsS = sqrt(mean(B.^2));
    rmsN = sqrt(mean((Bn-B).^2));
    
    title([num2str( mean(sqrt(B)) ),'-',...
           num2str( sqrt(mean(B)) ),'-',...
           num2str( rmsS/rmsN )])
    plot(Bn)
    plot(B)

end