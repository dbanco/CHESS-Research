% Data waterfall
figure
load("C:\Users\dpqb1\Documents\Dissertation\ms_sim1\output_j1_sig_1.00e-02_lam1_6.00e-02_lam2_0.00e+00.mat")
waterfall(flip(squeeze(outputs.y)',1))
xlabel('\eta','FontSize',18)
ylabel('time','FontSize',18)
zlabel('Intensity','FontSize',18)


%% Recon waterfall
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

figure
waterfall(flip(squeeze(Yhat)',1))
colorbar()

% xlabel('\eta','FontSize',18)
% ylabel('time','FontSize',18)
% zlabel('Intensity','FontSize',18)