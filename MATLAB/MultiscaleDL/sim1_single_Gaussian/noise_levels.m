sigmas = 0:0.01:0.1;
meanSNR = zeros(numel(sigmas),1);
for n = 2:numel(sigmas)
    [y,y_true,N,M,T] = gaus_example_multiscale_dl(sigmas(n));
    
    SNR = zeros(T,1);
    nPwr = zeros(T,1);
    sigPwr = zeros(T,1);
    for t = 1:T
        SNR(t) = norm(y_true(:,t))/norm(y(:,t)-y_true(:,t));
        nPwr(t) = norm(y(:,t)-y_true(:,t));
        sigPwr(t) = norm(y_true(:,t));
    end
    
    meanSNR(n) = mean(SNR);
end


NN = numel(sigmas);
trueErr = zeros(NN,1);
dataErr = zeros(NN,1);
noiseNorm = zeros(NN,1);
noiseNorm2 = zeros(NN,1);
outDir = "C:\Users\dpqb1\Documents\Outputs2024\";
lam_sel = [0.06,0.06,0.2,0.3,0.3,0.4,0.5,0.5,0.5,0.5,0.5];
for n = 2:NN
    gausDir = "gaus_example_8_8_24_X0_D0_V00_sig_"+num2str(n);
    dataFile = sprintf("output_j1_sig_%1.2s_lam1_%1.2s_lam2_0.00e+00",sigmas(n),lam_sel(n));
    [y,y_true,N,M,T] = gaus_example_multiscale_dl(sigmas(n));
    load(fullfile(outDir,gausDir,dataFile))
    D = outputs.D;
    N = outputs.N;
    M = outputs.M;
    y = outputs.y;
    X = outputs.X;
    scales = outputs.scales;
    center = (M+1)/2;
    
    AD = reSampleCustomArrayCenter(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);

    trueErr(n) = norm(y_true-Yhat);
    dataErr(n) = norm(squeeze(outputs.y)-Yhat);
    noiseNorm(n) = norm(randn(N,T)*sigmas(n));
    noiseNorm2(n) = norm(y_true-squeeze(outputs.y));
end

figure()
hold on
plot(meanSNR(2:NN),trueErr(2:NN),'o-')
plot(meanSNR(2:NN),dataErr(2:NN),'x-')
plot(meanSNR(2:NN),noiseNorm(2:NN),'s-')
plot(meanSNR(2:NN),noiseNorm2(2:NN),'x-')
xlabel('SNR','Fontsize',14)
ylabel('Error','Fontsize',14)
legend('$\|\hat{{\bf b}}-{\bf f}\|_2$','$\|\hat{{\bf b}}-{\bf b}\|_2$',...
    '$\|\mathcal{N}(0,s)\|_2$','$\|{\bf b}-{\bf f}\|_2$',...
    'interpreter','latex','Fontsize',14)

% close all
% figure(1)
% plot(SNR)
% figure(2)
% plot(nPwr)
% figure(3)
% plot(sigPwr)
% 
% norm(randn(N,1)*sigmas(2))

