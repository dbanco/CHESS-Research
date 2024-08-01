% Show simulated data

load("C:\Users\dpqb1\Documents\Dissertation\ms_sim2\output_j1_sig_1.00e-02_lam1_6.00e-02_lam2_0.00e+00.mat")
suffix = 'indep';
% load("C:\Users\dpqb1\Documents\Dissertation\ms_sim2\output_j10_sig_1.00e-02_lam1_6.00e-02_lam2_5.00e-01.mat")
% suffix = 'of';
figDir = "C:\Users\dpqb1\Documents\Dissertation\ms_sim2";

N = outputs.N;
M = outputs.M;
D = gather(outputs.D);  
suffix = 'true';
[y,y_true,K,J,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center(0);


D(:,:,2) = circshift(D(:,:,2),4);
norm(D(:)-Dtrue(:))/norm(Dtrue(:))

figure(1)
plot(D(:,:,1))
hold on
plot(Dtrue(:,:,1))

figure(21)
plot(D(:,:,2))
hold on
plot(Dtrue(:,:,2))

outputs.D = Dtrue;
X = gather(outputs.X);
y = gather(outputs.y);
scales = outputs.scales;
center = (M+1)/2;

generateFiguresToy1zpad_center(figDir,outputs,suffix,[4,4]);
