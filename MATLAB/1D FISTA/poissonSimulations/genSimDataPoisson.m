function [Bn,B,theta_stds,rel_err] = genSimDataPoisson(N,T,alpha,sim)
% Generate data
if strcmp(sim,'linear')
    theta_stds = linspace(1,15,T)';
elseif strcmp(sim,'anomaly')
    theta_stds = [7*ones(1,T/2),12*ones(1,T/2)]';
end

B = zeros(N,T);
Bn = zeros(N,T);
rel_err = zeros(T,1);
for t = 1:T
    close all
    b = gaussian_basis_1D( N, N/2, theta_stds(t)^2);
    rms = sqrt(sum(b.^2)/N);
    b = b*500/alpha;
    bn = poissrnd(b);
    B(:,t) = b/500*alpha/rms/3;
    Bn(:,t) = bn/500*alpha/rms/3;
    rel_err(t) = norm(b-bn)/norm(b);
end
end

