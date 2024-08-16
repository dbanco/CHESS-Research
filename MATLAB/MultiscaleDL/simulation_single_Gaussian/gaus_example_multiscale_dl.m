function [yn,y,N,M,T] = gaus_example_multiscale_dl(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 50;
N = 105; M = 105;

if nargin < 1
    sigma = 0.01;
end

% Position in time
position = 15+[linspace(20,50,20),30+20.*linspace(1,0.1,30).*cos((0:29)/3)];
width = (70-fliplr(position))/3;

% figure(1)
% plot(position)
% hold on
% plot(width)

for t = 1:T
    y(:,t) = gaussian_basis_wrap_1D(N,position(t),width(t),'2-norm');
end

yn = y + randn(N,T)*sigma;
% figure(2)
% imagesc(y)

end