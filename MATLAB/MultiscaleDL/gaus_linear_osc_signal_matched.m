function [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 150;
N = 300; M = 150;
if nargin < 1
    sigma = 0.015;
end

% Model Setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],8,100, 1/8);
scales{2} = genRationals([0;1],[1;1],8,100, 1/8);
Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
J = numel(scales);
KJ = sum(Uarray);
opt = [];
opt.DictFilterSizes = [1,1;...
                       M,M];

% Width in time
% sig = 2*linspace(2,10,T);
% width = 4*round(linspace(2,14,T));
pFactor = 3.5;
tt = (1:T)/(T/20)/pFactor;
sig = round(10.5*sin(tt) + 11.5);
width = round(10.5*cos(tt) + 11.5)+22;

sigMax = 18;
widthMax = 80;

% Define dictionary atoms
scaling = '2-norm';
Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = gaussian_basis_wrap_1D(M,80,sigMax,scaling);
Dtrue(1,21:widthMax+20,2) = linspace(0,10,widthMax)/norm(linspace(0,10,widthMax));
ADtrue = reSampleCustomArray(N,Dtrue,scales);
ADf = fft2(ADtrue);

% Position in time
minPos = 35+40;
maxPos = 107+50;
amp = (maxPos-minPos)*0.5;
Pos  = mod(round( amp*sin( 4*pi/T*(1:T)/pFactor + pi/2)          + amp + minPos)+300,N)+1;%-2*amp;%+4;
Pos2 =     round( amp*sin( 4*pi/T*(1:T)/pFactor + pi + pi/2) + amp + minPos)-30;%+2*amp;%-4;
Xtrue = zeros(1,N,KJ,T);
for t = 1:T
    Xtrue(1,Pos(t),sig(t),t) = 1;
    Xtrue(1,Pos2(t),width(t),t) = 1;
end
Xf = fft2(Xtrue);
y = ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric');
yn = y + randn(1,N,1,T)*sigma;
yn = reshape(yn,[1,N,T]);

% Reduce data to a time subset
% trange = 1:60;
% yn = yn(:,:,trange);
% y = y(:,:,trange);
% Xtrue = Xtrue(:,:,:,trange);
% T = numel(trange);

% imagesc(squeeze(yn))

end