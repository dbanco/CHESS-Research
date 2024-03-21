% Show simulated data

load("C:\Users\dpqb1\Documents\Dissertation\ms_sim2\output_j1_sig_1.00e-02_lam1_6.00e-02_lam2_0.00e+00.mat")

D = outputs.D;
N = outputs.N;
M = outputs.M;
y = squeeze(outputs.y)';
X = outputs.X;
scales = outputs.scales;
center = (M+1)/2;
AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
Yhat = gather(Yhat');
Yhat = Yhat+ 0.000000001*randn(size(Yhat));

azimuth = -17;
elevation = 45;
fig_pos = [600  100  500 300];
fig_pos2 = fig_pos+[0 400 0 0];

ff1 = figure(1);
waterfall(flip(y,1))
% xlabel('\eta','FontSize',18)
% ylabel('time','FontSize',18)
% zlabel('Intensity','FontSize',18)
xlim([0,55])
ylim([0,30])
zlim([-0.1,1.1])
ff1.Position = fig_pos;
view(azimuth, elevation);

ff2 = figure(2);
waterfall(flip(Yhat,1))
% xlabel('\eta','FontSize',18)
% ylabel('time','FontSize',18)
% zlabel('Intensity','FontSize',18)
xlim([0,55])
ylim([0,30])
zlim([-0.1,1.1])
ff2.Position = fig_pos2;
view(azimuth, elevation);


ff3 = figure(3);
waterfall(flip(Yhat,1))
% xlabel('\eta','FontSize',18)
% ylabel('time','FontSize',18)
% zlabel('Intensity','FontSize',18)
xlim([0,55])
ylim([0,30])
zlim([-0.1,1.1])
ff3.Position = fig_pos2;
view(azimuth, elevation);
colorbar()


load("C:\Users\dpqb1\Documents\Dissertation\ms_sim2\output_j10_sig_1.00e-02_lam1_6.00e-02_lam2_5.00e-01.mat")

D = outputs.D;
N = outputs.N;
M = outputs.M;
y = squeeze(outputs.y)';
X = outputs.X;
scales = outputs.scales;
center = (M+1)/2;
AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
Yhat = gather(Yhat');
Yhat = Yhat+ 0.000000001*randn(size(Yhat));


ff4 = figure(4);
waterfall(flip(Yhat,1))
% xlabel('\eta','FontSize',18)
% ylabel('time','FontSize',18)
% zlabel('Intensity','FontSize',18)
xlim([0,55])
ylim([0,30])
zlim([-0.1,1.1])
ff4.Position = fig_pos2;
view(azimuth, elevation);


simDir = 'C:\Users\dpqb1\Documents\Dissertation\ms_sim2\';
saveas(ff1,fullfile(simDir,'sim2_data.png'))
saveas(ff2,fullfile(simDir,'sim2_recon.png'))
saveas(ff3,fullfile(simDir,'sim2_colorbar.png'))
saveas(ff4,fullfile(simDir,'sim2_recon_of.png'))

