lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];

sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'sim2_gaussian_tooth_matched','sim2_gaussian_tooth_unmatched',...
            'sim2_gaussian_tooth_matched2','sim2_gaussian_tooth_unmatched2',...
            'dissertation','dissertation_long','dissertation_long_separate',...
            'pseudo-voigt_unmatched'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s0 = 8;
s1 = 2;
s2 = 1;
s3 = 2;
s4 = 1;
s5 = 1;
dataset = datasets{s0};
opt.Penalty = penalties{s1};
opt.coefInit = xinits{s2};
opt.dictInit = dinits{s3};
opt.Dfixed = dfixes{s4};
opt.Recenter = recenters{s5};
opt.Xfixed = 0;
topDir = ['E:\MCDLOF_processing\Outputs_7_24indep_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

% criterion = 'discrepancy';
% criterion = 'truth_errorr';
% criterion = 'relaxed discrepancy';
%criterion = 'discrepancy range';
criterion = 'discrepancy range2';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,1);
objectives = cell(NN,1);

useMin = 0;
relax_param = 0.004; % for discrep range2
sig_ind = 1:6;
for n = sig_ind
    inDir = [topDir,'\results_trial_1_sig_',num2str(n)];
    resultFile = [topDir,'\results_sig_',num2str(n),'.mat'];
    if exist(resultFile,'file')
        load(resultFile)
    else
        results = compute_metrics(inDir,sigmas(n),dataset,useMin,true);
        save(resultFile,'results')
    end
    [lambda_all,objective] = param_select_mcdl(results,criterion,sigmas(n),dataset,relax_param,fig_num);
    selected_lam_all_vec(n,:) = lambda_all;
    selected_inds(n,1) = find(selected_lam_all_vec(n,1) == lambdaVals);
    selected_inds(n,2) = find(selected_lam_all_vec(n,2) == lambdaOFVals);
    selected_inds(n,3) = find(selected_lam_all_vec(n,3) == lambdaHSVals);
    j_s = selected_inds(n,1);
    j_of = selected_inds(n,2);
    j_hs = selected_inds(n,3);
    if makeFigures
        outputs = loadOutputFile(inDir,selected_inds(n,:));
        suffix = sprintf('_j%i_%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e_lam3_%0.2e',...
                      j_s,j_of,j_hs,sigmas(n),outputs.lambda,outputs.lambda2,outputs.lambda3);
        psFigDir = fullfile(topDir,['ps_',criterion,'_',num2str(relax_param),...
                                          '_useMin',num2str(useMin)]);
        generateFiguresToy1zpad_center(psFigDir,outputs,suffix,[4,8],useMin);
        close all
    end
    objectives{n} = objective;
end
% LcurveFile = fullfile(topDir,'l-curve_plot.png'); 
% fig = gcf;
% fig.Position = [100 100 1400 800];
% saveas(gcf, LcurveFile);
% removeWhiteSpace(LcurveFile)

% %% Compute meanSNR
% [meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);
% 
% %% Compute errors
% dirStartS = 'results_1';
% [noiseNorm,trueErrS,dataErrS,~,~] = simError(y_true,sigmas,sig_ind,topDir,dirStartS,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals);
% 
% figure()
% hold on
% plot(meanSNR,noiseNorm,'s-')
% plot(meanSNR,trueErrS,'o-')
% xlabel('SNR','Fontsize',14)
% ylabel('Error','Fontsize',14)
% legend('$\|{\bf w}\|_2$','$\|\hat{{\bf b}}-{\bf f}\|_2$',...ub hbh h                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%     '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...
%     'interpreter','latex','Fontsize',14)
% 
% 
% %% Next copy figures associated with selected parameters to a folder
% pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim1_lam_s_6_26indep',...
%            '_relax_',num2str(relax_param),...
%            '_useMin1_',criterion,'_',opt.Penalty,'_',dirStartS,'.pptx'];
% titleStr = ['Sim 1 Param Select,',dirStartS];
% % createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)
% createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,sigmas,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind,objectives)