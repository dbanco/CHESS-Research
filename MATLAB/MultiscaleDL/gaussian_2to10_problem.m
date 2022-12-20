function [y,N,M,T] = gaussian_2to10_problem
%% Construct 1D test problem Box and Gaussian
T = 60;
N = 128; M = N;
K = 1;
U = 3;
scaling = '2-norm';

% Width in time
sig = linspace(2,10,T);

% Position in time
minPos = 28;
maxPos = 100;
amp = (maxPos-minPos)/2;
Pos = round( amp*sin( 2*pi*(1:T)/(T/2)) + amp + minPos);

y = zeros(1,N,T);
X_true = zeros(1,N,K*U,T);
for t = 1:T
    d1 = gaussian_basis_wrap_1D(N,1,sig(t),scaling);
    y(:,:,t) = y(:,:,t) + circshift([d1', zeros(1,N-M)],Pos(t));
end

end