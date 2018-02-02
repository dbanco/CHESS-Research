%% Load data
load('D:\CHESS_data\al7075_311_alpha_beta_part3\fista_fit_0_35_task_31.mat')
%% Compute and plot weighting functions
P.betap = 1;
P.alphap = 1;
weight = zeros(numel(P.var_theta),numel(P.var_rad));
for t = 1:numel(P.var_theta)
    for r = 1:numel(P.var_rad)
        weight(t,r) = P.alphap*2.^(-P.weight*P.betap*P.dtheta/sqrt(P.var_theta(t))*P.drad/sqrt(P.var_rad(r)));
    end
end

mesh(weight)
xlabel('\sigma^2_\theta')
ylabel('\sigma^2_r')