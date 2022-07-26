function [B,theta_stds] = genLinearPoisson(N,T,alpha)
% Generate data
theta_stds = linspace(1,15,T);
B = zeros(N,T);
for t = 1:T
    close all
    b = gaussian_basis_1D( N, N/2, theta_stds(t)^2);
    rms = sqrt(sum(b.^2)/N);
    b = b/rms*100;
    bn = poissrnd(b) + poissrnd(alpha*ones(N,1));
    B(:,t) = bn;
%         if t == 30
%             plot(b);hold on;plot(bn);
%         end
end
end

