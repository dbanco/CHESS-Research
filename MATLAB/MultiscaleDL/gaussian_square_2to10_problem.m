function [y,N,M,T] = gaussian_square_2to10_problem
%% Construct 1D test problem Box and Gaussian
T = 60;
N = 160+140; M = N;

scaling = '2-norm';

% Width in time
sig = linspace(2,10,T);
wid1 = 2*round(linspace(2,14,T));
wid2 = 2*round(linspace(1,7,T));

% Position in time
minPos = 35+70;
maxPos = 107+70;
amp = (maxPos-minPos)/2;
Pos  = round( amp*sin( 4*pi/T*(1:T))      + amp + minPos);
Pos2 = round( amp*sin( 4*pi/T*(1:T) + pi) + amp + minPos);

y = zeros(1,N,T);
for t = 1:T
    d1 = gaussian_basis_wrap_1D(N,1,sig(t),scaling);
    y(:,:,t) = y(:,:,t) + circshift([d1', zeros(1,N-M)],Pos(t));
end
for t = 1:T
    d2 = zeros(1,N);
    d2(1:wid1(t)) = 2;
    d2(1:wid2(t)) = 1;
    d2 = d2/norm(d2);
    y(:,:,t) = y(:,:,t) + circshift([d2, zeros(1,N-M)],Pos2(t));
end
imagesc(squeeze(y))
end