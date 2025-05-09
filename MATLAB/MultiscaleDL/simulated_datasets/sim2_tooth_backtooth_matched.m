function [yn,y,N,M,T,Xtrue,Dtrue] = sim2_tooth_backtooth_matched(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 30;
N = 55; 
M = 55;

if nargin < 1
    sigma = 0;
%     sigma = 0.015;
end

% Model Setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,8, 1/8);
scales{2} = genRationals([0;1],[1;1],16,8, 1/8);
J = size(scales{1},2);
KJ = K*J;
center = (M+1)/2;

tooth_width = 18;
Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = flip(tooth(M,tooth_width));
Dtrue(1,:,2) = tooth(M,tooth_width);

ADtrue = reSampleCustomArrayCenter(N,Dtrue,scales,center);
ADpad = padarray(ADtrue,[0 M-1 0 0],0,'post');
ADf = fft2(ADpad);

Pos = round(linspace(8,43,T)) + center;
Pos2 = round(linspace(43,8,T)) + center;

sig = round(linspace(1,8,T));
width = round(linspace(8,1,T));

Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,Pos(t),sig(t),t) = 1;
    Xtrue(1,Pos2(t),width(t)+J,t) = 1;
end

Xf = fft2(Xtrue);
y = squeeze(unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre'));

yn = y + randn(N,T)*sigma;

% figure(2)
% imagesc(yn)
% 
% figure(3)
% subplot(3,1,1)
% plot(yn(:,1),'-o')
% subplot(3,1,2)
% plot(yn(:,12),'-o')
% subplot(3,1,3)
% plot(yn(:,30),'-o')
end