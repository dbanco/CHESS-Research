position = [linspace(20,50,20),30+20.*linspace(1,0.1,30).*cos((0:29)/3)];
width = (85-fliplr(position))/3;
close all
figure(44)
plot(width)

width = width-min(width);
width = width*25/max(width);
width = width+1;

figure(44)
hold on
plot(width)

load("C:\Users\dpqb1\Documents\Dissertation\ms_sim1\output_j1_sig_1.00e-02_lam1_6.00e-02_lam2_0.00e+00.mat")

generateFiguresToy1zpad_center([],outputs,[],[4,8])
figure(3)
hold on
plot(width,'--','LineWidth',1,'Color','red')

% y = squeeze(outputs.y);
% figure(78)
% subplot(2,1,1)
% plot(y(:,5))
% subplot(2,1,2)
% plot(y(:,45))

% y_true = squeeze(outputs.y_true);
% figure(79)
% subplot(2,1,1)
% plot(y_true(:,5))
% subplot(2,1,2)
% plot(y_true(:,45))

%%

figure(80)
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
err = sum((squeeze(y)-Yhat).^2,1);
plot(err)


close all
figure(91)
atom = 0.92*circshift(AD(:,:,25),12);
hold on
plot(y(:,:,20),'LineWidth',2)
plot(Yhat(:,20),'LineWidth',2)
plot(atom(:,1:N),'LineWidth',2,'LineStyle',':')


%%
figure(911)
plot(AD(:,:,32),'LineWidth',2,'LineStyle',':')
