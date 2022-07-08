len = 2^11;
h = [4  -5  3  -4  5  -4.2   2.1   4.3  -3.1   5.1  -4.2];
t = [0.1  0.13  0.15  0.23  0.25  0.40  0.44  0.65  0.76  0.78  0.81];
h  = abs(h);
w  = 0.01*[0.5 0.5 0.6 1 1 3 1 1 0.5 0.8 0.5];
tt = linspace(0,1,len);
xref = zeros(1,len);
for j=1:11
    xref = xref+(h(j)./(1+((tt-t(j))/w(j)).^4));
end
% Add zero-mean white Gaussian noise with a variance of 0.25.
rng default
x = xref + 0.5*randn(size(xref));
plot(x)
axis tight

%% Generate example data
% N = 501;
% numSpots = 5;
% b = zeros(N,1);
% amplitude = rand(numSpots,1)*5;
% mean_param = rand(numSpots,1)*N;
% var_param = rand(numSpots,1)*N/10;
% for i = 1:numSpots
%    b = b + amplitude(i)*gaussian_basis_1D(N, mean_param(i), var_param(i));
% end
% 
% % Add noise
% b_n = b + 0.5*randn(N,1);
% xref = b
% x = b_n


figure(1)
origmode = dwtmode('status','nodisplay');
dwtmode('per','nodisplay')
xd = wdenoise(x,3,'Wavelet','sym4',...
    'DenoisingMethod','UniversalThreshold','NoiseEstimate','LevelIndependent');
plot(xd)
axis tight
hold on
plot(xref,'r')
legend('Denoised','Reference')


figure(2)
N = 8;
[C,L] = wavedec(x,N,'sym4')
n = 1;
D = detcoef(C,L,n)
subplot(2,1,1)
plot(C)
subplot(2,1,2)
plot(D)

figure(3)
[LO_D,HI_D,LO_R,HI_R] = wfilters('sym8')
subplot(2,2,1)
plot(LO_D)
title('Low decomp')
subplot(2,2,2)
plot(HI_D)
title('High decomp')
subplot(2,2,3)
plot(LO_R)
title('Low recon')
subplot(2,2,4)
plot(HI_R)
title('High recon')
