function [Bn,B,theta_stds,rel_err] = genAnomalyPoisson(N,T,alpha)
% Generate data
theta_stds = [7*ones(1,T/2),12*ones(1,T/2)];
B = zeros(N,T);
Bn = zeros(N,T);
rel_err = zeros(T,1);
for t = 1:T
    close all
    b = gaussian_basis_1D( N, N/2, theta_stds(t)^2);
    b = b*500/alpha;
    bn = poissrnd(b);
    rms = sqrt(sum(b.^2)/N);
    B(:,t) = b/rms;
    Bn(:,t) = bn/rms;
    rel_err(t) = norm(b-bn)/norm(b);
end
end

