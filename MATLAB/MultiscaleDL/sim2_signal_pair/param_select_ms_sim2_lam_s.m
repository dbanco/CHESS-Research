lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

tradeoff_s = 0.4;
tradeoff_of = 1;
scaleP = [0.4,5.4,9.48,116,0,100];

% criterion = 'discrepancy';
criterion = 'truth_error';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);

% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_10_Dtrue_Xtrue';
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_31_Dtrue1_Xzeros0';
% dirStartS = 'steps_matched_results';
fig_num = 22;

dataset = 'sim2_gaussian_tooth_matched';
testType = 'Drand0_Xzeros0';
topDir = ['E:\Outputs_sim2_lam_all2_',testType];
dirStartS = 'sim2_gaussian_tooth_matched_log_results';

sig_ind = 2:4;
for n = 2:4
    inDir = [topDir,'\',dirStartS,'_sig_',num2str(n)];
     % True error
    [~,y_true,~,~,~] = sim_switch_multiscale_dl(sigmas(n),dataset);
    [lambda_s_sel,j_s] = param_select_lambda_s(inDir,tradeoff_s,scaleP,22,criterion,sigmas(n),y_true);
    selected_lam_s_vec(n) = lambda_s_sel;
end
LcurveFile = fullfile(topDir,['l-curve_plot',dirStartS,'.png']); 
fig = gcf;
fig.Position = [100 100 1400 800];
saveas(gcf, LcurveFile);
removeWhiteSpace(LcurveFile)

%% Compute meanSNR
[meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);

%% Compute errors
[noiseNorm,trueErrS,dataErrS,~,~] = simError(y_true,sigmas,sig_ind,topDir,dirStartS,selected_lam_s_vec,lambdaVals);

figure()
hold on
plot(meanSNR,noiseNorm,'s-')
plot(meanSNR,trueErrS,'o-')
xlabel('SNR','Fontsize',14)
ylabel('Error','Fontsize',14)
legend('$\|{\bf w}\|_2$','$\|\hat{{\bf b}}-{\bf f}\|_2$',...ub hbh h                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...
    'interpreter','latex','Fontsize',14)


%% Next copy figures associated with selected parameters to a folder
pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim2_lam_s',criterion,'_',testType,'_',dirStartS,'.pptx'];
titleStr = ['Sim 2 Recovery, ',testType,', ',dirStartS];
createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)

