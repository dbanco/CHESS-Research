%% Testing natuer of 1-norm fit
num_theta = 161;
az_mean1 = 81;
az_std1 = 10;
truth = gaussian_basis_1D_norm2( num_theta, az_mean1, az_std1^2);

figure(1)
plot(truth)

lam = 0.1;
eps = logspace(-4,-1,500);
cost= zeros(size(eps));
for i = 1:numel(eps)  
    fit = gaussian_basis_1D_norm2( num_theta, az_mean1, (az_std1 + eps(i)).^2);
    x = (1-eps);
    cost(i) = norm(truth-fit*x) + lam*x;
end

figure(2)
plot(cost)