function [alpha_val,theta_stds] = genPoissonScales(N,T,levels,sim)
%% Comptue error for different data scales
alphas = linspace(1,500,1000);
A = numel(alphas);
rel_errs = zeros(A,1);

if strcmp(sim,'linear')
    for a = 1:A
        [~,~,theta_stds,rel_err] = genLinearPoisson(N,T,alphas(a));
        rel_errs(a) = mean(rel_err);
    end
elseif strcmp(sim,'anomaly')
    for a = 1:A
        [~,~,theta_stds,rel_err] = genAnomalyPoisson(N,T,alphas(a));
        rel_errs(a) = mean(rel_err);
    end
end
% plot(alphas,rel_errs)

%% Choose scales to get particular error
alpha_val = zeros(numel(levels),1);
for i = 1:numel(levels)
    [errDist,ind] = min(abs(rel_errs-levels(i)));
    alpha_val(i) = alphas(ind);
    err_val(i) = rel_errs(ind);
end

%% Show selected scales/err levels
% for a = 1:numel(levels)
%     [Bn,B,~,~] = genLinearPoisson(N,T,alpha_val(a));
%     
%     t=30;
%     figure(1)
%     hold on
%     plot(B(:,t))
%     plot(Bn(:,t))
%     title(num2str(alpha_val(a)))
%     pause()
% end