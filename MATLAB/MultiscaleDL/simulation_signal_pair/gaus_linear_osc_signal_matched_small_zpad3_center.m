function [yn,y,K,J,N,M,T,Xtrue,Dtrue,scales] = gaus_linear_osc_signal_matched_small_zpad3_center(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 30;
N = 55; M = 45;
center = (M+1)/2;
if nargin < 1
    sigma = 0.02;
%     sigma = 0.015;
end

% Model Setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,8, 1/6);
scales{2} = genRationals([0;1],[1;1],8,8, 1/6);
J = size(scales{1},2);
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
tt = ((1:T))/5;
sig = round(   (J-1)/2*sin(tt) + (J-1)/2+1 );
width = round( (J-1)/2*cos(tt) + (J-1)/2+1 )+ J;
 
sigMax = 4;
widthMax = 20;
 
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
amp = 15;
Pos  = M-1-1 + mod(round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi/2 )      + amp + minPos1),N);
Pos2 =     M-1-1 + round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi + pi/2 ) + amp + minPos2);
Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,Pos(t),sig(t),t) = 1;
    Xtrue(1,Pos2(t),width(t),t) = 1;
end
Xf = fft2(Xtrue);
y = unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre');
yn = y + randn(1,N,1,T)*sigma;
yn = reshape(yn,[1,N,T]);

% snr = norm(y(:))/norm(y(:)-yn(:))
% Reduce data to a time subset
% trange = 1:60;
% yn = yn(:,:,trange);
% y = y(:,:,trange);
% Xtrue = Xtrue(:,:,:,trange);
% T = numel(trange);

figure
imagesc(squeeze(yn))

f1 = figure;
hold on
for i = 1:J
%     subplot(7,7,i)
    plot(real(ADtrue(:,:,i)),'Linewidth',1)
    set(gca, 'XtickLabel','')
    set(gca, 'FontSize', 16)
end
f1.Position = [1 100 900 500];

f2 = figure;
hold on
for i = J+1:J+J
%     subplot(7,7,i)
    plot(real(ADtrue(:,:,i)),'Linewidth',1)
    set(gca, 'XtickLabel','')
    set(gca, 'FontSize', 16)
end
f2.Position = [1 100 900 500];

figure
imagesc(squeeze(sum(Xtrue,[1,2])))

end