function [yn,y,N,M,T,Xtrue,Dtrue] = peak_example_steps_matched_multiscale_dl(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 64;
N = 105; M = 105;

if nargin < 1
    sigma = 0.01;
end

width = ones(T,1);
val = 0;
for t = 1:T
    if mod(t-1,4) == 0
        val = val + 1;
    end
    width(t) = val;

end

% Model Setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
KJ = K*J;
center = (M+1)/2;

load('D:\MMPAD_data_nr1\ring4_zero\mmpad_img_50.mat')
peakData = sum(polar_image(:,185:229),1);
peakData = peakData-min(peakData);
Dtrue = zeros(1,M,K);
Dtrue(1,31:75,1) = peakData;
ADtrue = reSampleCustomArrayCenter(N,Dtrue,scales,center);
ADpad = padarray(ADtrue,[0 M-1 0 0],0,'post');
ADf = fft2(ADpad);

Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,M,width(t),t) = 1;
end

Xf = fft2(Xtrue);

% figure(1)
% sigma_true = zeros(J,1);
% for j = 1:J
%     mdl = fittype('gauss1');
%     f = fit((1:N)',ADtrue(1,:,j)',mdl);
%     sigma_true(j) = f.c1;
%     ADtrue(1,:,j) = f(1:N)/norm(f(1:N));
%     subplot(2,J/2,j)
%     plot(f,(1:N)',ADtrue(1,:,j))
%     if j ~= J
%         lgd = findobj('type', 'legend');
%         delete(lgd)
%     end
% end

y = squeeze(unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre'));

% y = zeros(N,T);
% 
% for t = 1:T/2
%     if t <= J
%         y(:,2*t) = gaussian_basis_wrap_1D(N,center,sigma_true(t),'2-norm');
%         y(:,2*t-1) = gaussian_basis_wrap_1D(N,center,sigma_true(t),'2-norm');
%         Xtrue(1,center,t,2*t) = 1;
%         Xtrue(1,center,t,2*t-1) = 1;
%     else
%         tt = mod(t,J+1)+1;
%         try
%             sig_half = (sigma_true(tt) + sigma_true(tt+1))/2;
%         catch
%             sig_half = sigma_true(tt) + 1;
%         end
% 
%         y(:,2*t) = gaussian_basis_wrap_1D(N,center,sig_half,'2-norm');
%         y(:,2*t-1) = gaussian_basis_wrap_1D(N,center,sig_half,'2-norm');
%         Xtrue(1,center,t-J,2*t) = 1;
%         Xtrue(1,center,t-J,2*t-1) = 1;
% 
%     end
% 
% end

yn = y + randn(N,T)*sigma;

figure(22)
imagesc(yn)

figure(33)
imagesc(y)

end