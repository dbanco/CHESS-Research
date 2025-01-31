lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

% criterion = 'discrepancy';
criterion = 'truth_error';

% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_10_Dtrue_Xtrue';
% topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_31_Dtrue1_Xzeros0';
% dirStartS = 'steps_matched_results';
fig_num = 22;

dataset = 'sim2_gaussian_tooth_unmatched2';
testType = 'Dflat0_Xzeros0';
topDir = ['E:\MCDLOF_processing\Outputs_',dataset,'_',testType];
dirStartS = [dataset,'_log_results'];

selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,1);
sig_ind = 2:4;
for n = sig_ind
    inDir = [topDir,'\',dirStartS,'_sig_',num2str(n)];
     % True error
    [~,y_true,~,~,~] = sim_switch_multiscale_dl(sigmas(n),dataset);
    [lambda_all,selInd,objective] = param_select_3D(inDir,fig_num,criterion,sigmas(n),y_true);
    selected_lam_all_vec(n,:) = lambda_all;
    selected_inds(n) = selInd;
end
LcurveFile = fullfile(topDir,['l-curve_plot',dirStartS,'.png']); 
fig = gcf;
fig.Position = [100 100 1400 800];
saveas(gcf, LcurveFile);
removeWhiteSpace(LcurveFile)

[meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);
true_error = objective.true_error(selected_inds(sig_ind));

figure()
hold on
plot(meanSNR,true_error,'s-')
plot(meanSNR,noiseError,'o-')
xlabel('SNR','Fontsize',14)
ylabel('Error','Fontsize',14)
legend('$\|{\bf w}\|_2$','$\|\hat{{\bf b}}-{\bf f}\|_2$',...ub hbh h                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...
    'interpreter','latex','Fontsize',14)


%% Next copy figures associated with selected parameters to a folder
pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim2_lam_s',criterion,'_',testType,'_',dirStartS,'.pptx'];
titleStr = ['Sim 2 Recovery, ',testType,', ',dirStartS];
createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,dirStartS,sigmas,...
    selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind)

