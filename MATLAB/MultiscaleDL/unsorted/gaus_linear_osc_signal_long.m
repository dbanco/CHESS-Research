function [yn,y,N,M,T] = gaus_linear_osc_signal_long(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 60;
N = 700; M = N;
if nargin < 1
    sigma = 0.015;
end
scaling = '2-norm';

% Width in time
% sig = 2*linspace(2,10,T);
% width = 4*round(linspace(2,14,T));
tt = (1:T)/3;
sig = round(8*sin(tt) + 16);
width = 3*round(12*sin(tt) + 16);

% Position in time
minPos = 35+90;
maxPos = 107+100;
amp = (maxPos-minPos);
Pos  = round( amp*sin( 4*pi/T*(1:T))      + amp + minPos)-50;%-2*amp;%+4;
Pos2 = round( amp*sin( 4*pi/T*(1:T) + pi) + amp + minPos)-50;%+2*amp;%-4;

y = zeros(N,T);
yn = zeros(N,T);
for t = 1:T
    d1 = gaussian_basis_wrap_1D(N,1,sig(t),scaling);
    y(:,t) = y(:,t) + circshift([d1; zeros(N-M,1)],Pos(t));
end
for t = 1:T
    d2 = zeros(N,1);
    d2(1:width(t)) = linspace(0,10,width(t));
%     d2(width(t)+1:2*width(t)-1) = d2((width(t)-1):-1:1);
    d2 = d2/norm(d2);
    y(:,t) = y(:,t) + circshift([d2; zeros(N-M,1)],Pos2(t));
end

% Add a second pair of signals
minPos = 35+90 + 350;
maxPos = 107+100 +350;
amp = (maxPos-minPos);
Pos  = round( amp*sin( 0.9*4*pi/T*(1:T) + pi/2)      + amp + minPos)-50;%-2*amp;%+4;
Pos2 = round( amp*sin( 0.8*4*pi/T*(1:T) + pi + pi/2) + amp + minPos)-50;%+2*amp;%-4;

for t = 1:T
    d1 = gaussian_basis_wrap_1D(N,1,sig(t),scaling);
    y(:,t) = y(:,t) + circshift([d1; zeros(N-M,1)],Pos(t));
end
for t = 1:T
    d2 = zeros(N,1);
    d2(1:width(t)) = linspace(0,10,width(t));
%     d2(width(t)+1:2*width(t)-1) = d2((width(t)-1):-1:1);
    d2 = d2/norm(d2);
    y(:,t) = y(:,t) + circshift([d2; zeros(N-M,1)],Pos2(t));
%     rms = sqrt(sum(y(:,t).^2)/N);
    bn = y(:,t) + randn(N,1)*sigma;
    yn(:,t) = bn;
end

% snr = norm( y(:) )/norm( y(:)-yn(:) )
imagesc(squeeze(yn))
yn = reshape(yn,[1,N,T]);
end