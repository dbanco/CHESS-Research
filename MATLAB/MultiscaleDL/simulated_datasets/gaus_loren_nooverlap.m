function [yn,y,N,M,T,Xtrue,Dtrue] = gaus_loren_nooverlap(sigma,plotFlag)
%% Construct 1D test problem Gaussian and linear 
T = 30;
N = 75; M = 31;
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

sig1 = round(   4*sin(0.1143*(1:T) + 2.2857) + 5.2);
sig2 = round( 3.8*sin(0.1143*(1:T) + 2.2857 + pi/2) + 5.2 )+ J;
 
sigMax = 4;
sig2Max = 4;
 
% Define dictionary atoms
Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = voigt_basis_wrap_1D(M,center,sigMax,0,'2-norm');
Dtrue(1,:,2) = voigt_basis_wrap_1D(M,center,sig2Max,1,'2-norm');
ADtrue = padarray(reSampleCustomArrayCenter3(N,Dtrue,scales,center),[0,M-1,0,0],0,'post');
ADf = fft2(ADtrue);

% Position in time
Pos =  round(15.9368*sin(0.0659*(1:T) + 1 ) + 18);
Pos2 = round(15.9368*sin(0.0659*(1:T) ) + 55);

Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,Pos(t),sig1(t),t) = 1;
    Xtrue(1,Pos2(t),sig2(t),t) = 1;
end
Xf = fft2(Xtrue);
y = unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre');
yn = y + randn(1,N,1,T)*sigma;
yn = reshape(yn,[N,T]);
y = reshape(y,[N,T]);

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
            i = i + 1;
        end
    end

    figure(6)
    window = 1:N+M-1;
    opt.HSiters = 100;
    plotOpticalFlow_sumT(Xtrue,K,opt,window)
end

end