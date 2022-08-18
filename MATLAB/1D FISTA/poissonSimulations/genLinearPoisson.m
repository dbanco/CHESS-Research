 function [Bn,B,theta_stds,rel_err] = genLinearPoisson(N,T,alpha)
% Generate data
theta_stds = linspace(1,15,T);
B = zeros(N,T);
Bn = zeros(N,T);
rel_err = zeros(T,1);
for t = 1:T
    close all
    b = gaussian_basis_1D( N, N/2, theta_stds(t)^2);
    rms = sqrt(sum(b.^2)/N);
    b = b*500/alpha;
    bn = poissrnd(b);
    B(:,t) = b;
    Bn(:,t) = bn;
    rel_err(t) = norm(b-bn)/norm(b);
end
end

