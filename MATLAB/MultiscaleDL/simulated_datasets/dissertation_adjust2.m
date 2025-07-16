function [yn,y,N,M,T,Xtrue,Dtrue] = dissertation_adjust2(sigma,plotFlag)
%% Construct 1D test problem Gaussian and linear 
T = 30;
N = 55; M = 31;
center = (M+1)/2;
if nargin < 1
    sigma = 0.01;
%     sigma = 0.015;
end
if nargin < 2
    plotFlag = false;
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
% pFactor = 3.5;
% tt = ((1:T)+20)/8.75;
% sigOLD =   round( 3.5*sin(tt) + 4.5 );
% widthOLD = round( 3.5*cos(tt) + 4.5 )+ J;3

sig = round(   4*sin(0.1143*(1:T) + 2.2857) + 5.2);
width = round( 3.8*sin(0.1143*(1:T) + 2.2857 + pi/2) + 5.2 )+ J;
 
sigMax = 4;
widthMax = 25;
 
% Define dictionary atoms
scaling = '2-norm';
Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = gaussian_basis_wrap_1D(M,center,sigMax,scaling);
i1 = round(center) - widthMax/2;
i2 = i1 + widthMax-1;
Dtrue(1,i1:i2,2) = linspace(0,10,widthMax)/norm(linspace(0,10,widthMax));
% Dtrue(1,1:widthMax,2) = linspace(0,10,widthMax)/norm(linspace(0,10,widthMax));
ADtrue = padarray(reSampleCustomArrayCenter3(N,Dtrue,scales,center),[0,M-1,0,0],0,'post');
ADf = fft2(ADtrue);

% Position in time
% minPos1 = 0;
% minPos2 = -5;
% amp = 15;
% PosOLD  = M-8 + mod(round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi/2 )  + amp + minPos1),N);
% Pos2OLD = M-8 + round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi + pi/2 ) + amp + minPos2);
Pos =  round(15.9368*sin(0.0659*(1:T) + 2.2947) + 48.1417);
Pos2 = round(15.9368*sin(0.0659*(1:T) - 0.8469) + 42.8583);

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

% snr = norm(y(:))/norm(y(:)-yn(:))
% % Reduce data to a time subset
% trange = 1:60;
% yn = yn(:,:,trange);
% y = y(:,:,trange);
% Xtrue = Xtrue(:,:,:,trange);
% T = numel(trange);
if plotFlag
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

    figure(6)
    window = 25:65;
    opt.HSiters = 100;
    plotOpticalFlow_sumT(Xtrue,K,opt,window)
end

end