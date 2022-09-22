function [Bn,B,theta_stds,rel_err,snr] = genSimDataPoisson2(N,T,alpha,sim)
% Generate data
if strcmp(sim,'linear')
    theta_stds = linspace(3,15,T)';
elseif strcmp(sim,'anomaly')
    theta_stds = [7*ones(1,T/2),12*ones(1,T/2)]';
end

B = zeros(N,T);
Bn = zeros(N,T);
rel_err = zeros(T,1);
snr = zeros(T,1);
for t = 1:T
    b = gaussian_basis_1D( N, N/2, theta_stds(t)^2);
    rms = sqrt(mean(b.^2));
    b = b*alpha;
    bn = poissrnd(b);
    B(:,t) = b/alpha;
    Bn(:,t) = bn/alpha;
    rel_err(t) = norm(b-bn)/norm(b);
    snr(t) = sqrt(mean(b.^2))/sqrt(mean((bn-b).^2));
end
end

