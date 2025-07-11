function [yn,y,N,M,T,Xtrue,Dtrue] = gaussian_tooth_matched_long2(sigma,plotFlag)
%% Construct 1D test problem Gaussian and linear 
T = 100;
N = 85; M = 35;
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
tt = linspace(0,3*pi,T);
sig = round(   -(J-1)/2*cos(tt) + (J-1)/2+1 );
width = round( (J-1)/2*cos(tt) + (J-1)/2+1 )+ J;
 
voigt_max_width = 5;
tooth_max_width = 26;

mixing = 1; % 1 = Gaussian, 0 = Lorentzian
 
% Define dictionary atoms
scaling = '2-norm';
Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = voigt_basis_wrap_1D(M,center,voigt_max_width,mixing,scaling);
i1 = round(center) - tooth_max_width/2;
i2 = i1 + tooth_max_width-1;
Dtrue(1,i1:i2,2) = linspace(0,10,tooth_max_width)/norm(linspace(0,10,tooth_max_width));
ADtrue = padarray(reSampleCustomArrayCenter3(N,Dtrue,scales,center),[0,M-1,0,0],0,'post');
ADf = fft2(ADtrue);

% Position in time
% minPos1 = 0;
% minPos2 = -5;
% amp = 15;
% Pos  = M-1-5 + mod(round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi/2 )      + amp + minPos1),N);
% Pos2 =     M-1-5 + round( amp*sin( 4*pi/50*((1:T)+9)/pFactor + pi + pi/2 ) + amp + minPos2);


amplitude = 20;

f1 = cos(linspace(0,2*pi,T));
f2 = -f1;

Pos1 = round(amplitude*f1 + M -center + (N-1)/2 );
Pos2 = round(amplitude*f2 + M-center + (N-1)/2 );

Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,Pos1(t),sig(t),t) = 1;
    Xtrue(1,Pos2(t),width(t),t) = 1;
end
Xf = fft2(Xtrue);
y = unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre');
yn = y + randn(1,N,1,T)*sigma;
yn = reshape(yn,[N,T]);
y = reshape(y,[N,T]);

snr = norm(y(:))/norm(y(:)-yn(:));

if plotFlag
    figure(2)
    imagesc(yn)
    
    figure(3)
    subplot(3,1,1)
    plot(yn(:,1),'-o')
    subplot(3,1,2)
    plot(yn(:,12),'-o')
    subplot(3,1,3)
    plot(yn(:,T),'-o')
    
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


end