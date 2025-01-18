lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

tradeoff_s = 0.4;
tradeoff_of = 1;
scaleP = [0.4,5.4,9.48,116,0,100];

criterion = 'discrepancy';
% criterion = 'truth_error';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);

% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_10_Dtrue_Xtrue';
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_31_Dtrue1_Xzeros0';
% dirStartS = 'steps_matched_results';
fig_num = 22;

dataset = 'steps_matched';
testType = 'Dflat0_Xzeros0';
topDir = ['E:\Outputs2024_12_3_',testType];
dirStartS = 'steps_matched_log_results';

for n = 1:NN
    inDir = [topDir,'\',dirStartS,'_sig_',num2str(n)];
     % True error
    [~,y_true,~,~,~] = gaus_example_switch_multiscale_dl(sigmas(n),dataset);
    [lambda_s_sel,j_s] = param_select_lambda_s(inDir,tradeoff_s,scaleP,22,criterion,sigmas(n),y_true);
    selected_lam_s_vec(n) = lambda_s_sel;
end
LcurveFile = fullfile(topDir,['l-curve_plot',dirStartS,'.png']); 
fig = gcf;
fig.Position = [100 100 1400 800];
saveas(gcf, LcurveFile);
removeWhiteSpace(LcurveFile)

%% Next plot the true recovery as a function of SNR
sigmas = 0:0.01:0.1;
meanSNR = zeros(numel(sigmas),1);
for n = 1:NN
    [y,y_true,N,M,T] = gaus_example_switch_multiscale_dl(sigmas(n),dataset);
    SNR = zeros(T,1);
    nPwr = zeros(T,1);
    sigPwr = zeros(T,1);
    for t = 1:T
        SNR(t) = norm(y_true(:,t))/norm(y(:,t)-y_true(:,t));
        nPwr(t) = norm(y(:,t)-y_true(:,t));
        sigPwr(t) = norm(y_true(:,t));
    end
    
    meanSNR(n) = mean(SNR);
end
%% Need to write a function to compute the various errors for all noise levels
% and that maybe just takes the file name as input so that it works on both
[noiseNorm,trueErrS,dataErrS,~,~] = simError(y_true,sigmas,topDir,dirStartS,selected_lam_s_vec,lambdaVals);

%     [noiseNorm,trueErrS,dataErrS,trueErrOF,dataErrOF] = simError(y_true,sigmas,outDirOF,dirStartS,...
%                                 selected_lam_s_vec,dirStartOF,selected_lam_of_vec,lambdaOFVals);
    
%%
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
% pptFile = 'C:\Users\dpqb1\Documents\MCDL Paper\recons_dicts_sim1_hs_2.pptx';
% titleStr = 'Sim 1 Recovery';
% createPowerpointSim(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,lambdaVals,dirStartOF,...
%     lambdaOFVals,lambdaHSVals,selected_lam_of_vec,selected_lam_s_vec)

pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\recons_1_17_',criterion,'_',testType,'_',dirStartS,'.pptx'];
titleStr = ['Sim 1 Recovery, ',testType,', ',dirStartS];
createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion)

