function [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps_multiscale_dl(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 50;
N = 105; M = 105;

if nargin < 1
    sigma = 0.01;
end

% Position in time
position = 15+[linspace(20,50,20),30+20.*linspace(1,0.1,30).*cos((0:29)/3)];
width = ones(T,1);
val = 2;
for t = 1:T
    width(t) = val;
    if mod(t,5) == 0
        val = val + 1;
    end
end

center = (M+1)/2;

position = round(position)+center;
width = round(width);

y = zeros(N,T);

for t = 1:T
    y(:,t) = gaussian_basis_wrap_1D(N,position(t)+center,width(t),'2-norm');
end

yn = y + randn(N,T)*sigma;

figure(22)
imagesc(yn)

end