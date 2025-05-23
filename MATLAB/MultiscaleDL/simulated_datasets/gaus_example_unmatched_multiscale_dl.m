function [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_unmatched_multiscale_dl(sigma,plotFlag)
%% Construct 1D test problem Gaussian and linear 
T = 50;
N = 105; M = 105;

if nargin < 1
    sigma = 0.01;
end
if nargin < 2
    plotFlag = false;
end

% Position in time
position = 15+[linspace(20,50,20),30+20.*linspace(1,0.1,30).*cos((0:29)/3)];
width = (70-fliplr(position))/3;

max_width = 10;

% Model Setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
KJ = K*J;
center = (M+1)/2;

Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = gaussian_basis_wrap_1D(M,center,max_width,'2-norm');
ADtrue = reSampleCustomArrayCenter(N,Dtrue,scales,center);

position = round(position)+center;
width = round(width);

% figure(1)
sigma_true = zeros(J,1);
for j = 1:J
    mdl = fittype('gauss1');
    f = fit((1:N)',ADtrue(1,:,j)',mdl);
    sigma_true(j) = f.c1;
    ADtrue(1,:,j) = f(1:N)/norm(f(1:N));
    % subplot(2,J/2,j)
    % plot(f,(1:N)',ADtrue(1,:,j))
    % if j ~= J
    %     lgd = findobj('type', 'legend');
    %     delete(lgd)
    % end
end

Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,position(t),width(t),t) = 1;
end
Xf = fft2(Xtrue);

ADpad = padarray(ADtrue,[0 M-1 0 0],0,'post');
ADf = fft2(ADpad);

y = squeeze(unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre'));
yn = y + randn(N,T)*sigma;

if plotFlag
    figure(2)
    imagesc(yn)
    
    figure(3)
    subplot(3,1,1)
    plot(yn(:,1),'-o')
    subplot(3,1,2)
    plot(yn(:,25),'-o')
    subplot(3,1,3)
    plot(yn(:,50),'-o')
    
    figure(4)
    vdf = squeeze(sum(squeeze(Xtrue),1));
    imagesc(vdf)
    
    figure(5)
    i = 1;
    for k = 1:2
        for j = 1:8
            subplot(2,8,i)
            plot(ADtrue(1,:,i))
            axis([0,N,0,0.9])
            i = i + 1;
        end
    end
end

end