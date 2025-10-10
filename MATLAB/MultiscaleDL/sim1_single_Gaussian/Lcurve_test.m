% Analyze regularization grid

% Load SNR = 12 for sim 2
load("E:\MCDLOF_processing\Outputs_10_6_reg_softmin_LBFGS_dissertation_adjust2_log_Dflat0_Xzeros0\results_sig_3.mat")

x = results.error;
y = results.reg_penalty;
z = results.log_penalty;

figure
plot3(x,y,z,'o')
xlabel('Error')
ylabel('Softmin Penalty')
zlabel('Log Penalty')

hold on

a1 = min((x));
a2 = min((y));
a3 = min((z));
% Compute origin distances
dist1 = log10(x - a1);
dist2 = log10(y - a2);
dist3 = log10(z - a3);
dist4 = dist1.^2 + dist2.^2 + dist3.^2;
[~,selInd] = min(dist4);

plot3(x(selInd),y(selInd),z(selInd),'s')

%% Plot lambda_reg vs reg_penalty
figure
loglog(results.lambda_vec(:,2),y,'o')
xlabel('Lambda2')
ylabel('Softmin Penalty')

% Plot error vs reg_penalty
figure
loglog(x,y,'o')
xlabel('Error')
ylabel('Softmin Penalty')
hold on

a1 = min((x));
a2 = min((y));

b1 = max((x));
b2 = max((y));

loglog(a1,a2,'o')

% Compute origin distances
dist1 = log10(x - a1);
dist2 = log10(y - a2);

dist3 = dist1.^2 + dist2.^2;

[~,selInd] = min(dist3);

loglog(x(selInd),y(selInd),'s')




