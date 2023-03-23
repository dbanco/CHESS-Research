function [yn,y,N1,N2,M,T] = gaus_linear_osc_signal_2D(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 200;
N1 = 60; N2= 350; M = N2;
if nargin < 1
    sigma = 0.015;
end
scaling = '2-norm';

% Width in time
% sig = 2*linspace(2,10,T);
% width = 4*round(linspace(2,14,T));
tt = (1:T)/(T/20);
sig = round(8*sin(tt) + 16);
width = 3*round(12*sin(tt) + 16);

% Position in time
minPos = 35+90;
maxPos = 107+100;
amp = (maxPos-minPos);
Pos  = round( amp*sin( 4*pi/T*(1:T))      + amp + minPos)-50;%-2*amp;%+4;
Pos2 = round( amp*sin( 4*pi/T*(1:T) + pi) + amp + minPos)-50;%+2*amp;%-4;

y = zeros(N1,N2,T);
yn = zeros(N1,N2,T);
for t = 1:T
    d1 = gaussian_basis_wrap_2D(N1,30,5,N2,1,sig(t),scaling);
    y(:,:,t) = y(:,:,t) + circshift([d1; zeros(N2-M,1)],[1 Pos(t)]);
end
for t = 1:T
    d2 = zeros(N1,N2,1);
    d2(1:15,1:width(t)) = repmat(linspace(0,10,width(t)), [15,1]);
%     d2(width(t)+1:2*width(t)-1) = d2((width(t)-1):-1:1);
    d2 = d2/norm(d2(:));
    y(:,:,t) = y(:,:,t) + circshift([d2; zeros(N2-M,1)],[23 Pos2(t)]);
%     rms = sqrt(sum(y(:,t).^2)/N);
    bn = y(:,:,t) + randn(N1,N2)*sigma;
    yn(:,:,t) = bn;
end
% snr = norm( y(:) )/norm( y(:)-yn(:) )
imagesc(squeeze(yn(:,:,20)))
yn = reshape(yn,[N1,N2,T]);
end