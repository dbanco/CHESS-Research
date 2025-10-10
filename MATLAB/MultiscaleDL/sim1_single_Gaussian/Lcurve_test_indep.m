% Analyze regularization grid

% Load SNR = 12 for sim 2
load("E:\MCDLOF_processing\Outputs_10_3_softmin_LBFGS_dissertation_adjust2_log_Dflat0_Xzeros0\results_sig_3.mat")

x = results.error;
y = results.log_penalty;

% Plot lambda_reg vs y
figure
loglog(results.lambda_vec(:,2),y,'o')
xlabel('Lambda2')
ylabel('Log Penalty')

% Plot error vs y
figure
loglog(x,y,'o')
xlabel('Error')
ylabel('Log Penalty')
hold on

b1 = max((x));
b2 = max((y));

% a1 = min((x));
% a2 = min((y));
val = 0.95;
a1 = min(x(x<b1*val));
a2 = min(y(x<b1*val));

loglog(a1,a2,'o')

% Compute origin distances
dist1 = log10(x - a1);
dist2 = log10(y - a2);

dist3 = dist1.^2 + dist2.^2;

[~,selInd] = min(dist3(x<b1*val));

loglog(x(selInd),y(selInd),'s')




