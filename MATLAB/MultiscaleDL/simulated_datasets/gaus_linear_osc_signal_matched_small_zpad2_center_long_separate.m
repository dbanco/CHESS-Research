function [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center_long_separate(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 100;
N = 55; M = 45;
center = (M+1)/2;
if nargin < 1
    sigma = 0.01;
%     sigma = 0.015;
end

% Model Setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,8, 1/6);
scales{2} = genRationals([0;1],[1;1],8,8, 1/6);
Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
J = size(scales{1},2);
KJ = sum(Uarray);
opt = [];
opt.DictFilterSizes = [1,1;...
                       M,M];

% Width in time
% sig = 2*linspace(2,10,T);
% width = 4*round(linspace(2,14,T));
pFactor = 3.5;
tt = ((1:T)+20)/(50/20)/pFactor;
sig = round(   (J-1)/2*sin(tt) + (J-1)/2+1 );
width = round( (J-1)/2*cos(tt) + (J-1)/2+1 )+ J;
 
sigMax = 4;
widthMax = 30;
 
% Define dictionary atoms
scaling = '2-norm';
Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = gaussian_basis_wrap_1D(M,center,sigMax,scaling);
Dtrue(1,6:widthMax+5,2) = linspace(0,10,widthMax)/norm(linspace(0,10,widthMax));
ADtrue = padarray(reSampleCustomArrayCenter(N,Dtrue,scales,center),[0,M-1,0,0],0,'post');
ADf = fft2(ADtrue);

% Position in time
minPos1 = 0;
minPos2 = -5;
amp = 1;
Pos  = M-11 + mod(round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi/2 )      + amp + minPos1),N);
Pos2 =     M-1+20 + round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi + pi/2 ) + amp + minPos2);
Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,Pos(t),sig(t),t) = 1;
    Xtrue(1,Pos2(t),width(t),t) = 1;
end
Xf = fft2(Xtrue);
y = unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre');
yn = y + randn(1,N,1,T)*sigma;
yn = reshape(yn,[N,T]);
y = reshape(y,[N,T]);

snr = norm(y(:))/norm(y(:)-yn(:))
% Reduce data to a time subset
% trange = 1:60;
% yn = yn(:,:,trange);
% y = y(:,:,trange);
% Xtrue = Xtrue(:,:,:,trange);
% T = numel(trange);

figure(2)
imagesc(yn)

figure(3)
subplot(3,1,1)
plot(yn(:,1),'-o')
subplot(3,1,2)
plot(yn(:,12),'-o')
subplot(3,1,3)
plot(yn(:,30),'-o')

figure(4)
vdf = squeeze(sum(squeeze(Xtrue),1));
imagesc(vdf)

figure(5)
i = 1;
for k = 1:2
    for j = 1:8
        subplot(2,8,i)
        plot(ADtrue(1,:,i))
        axis([0,60,0,0.9])
        i = i + 1;
    end
end

end