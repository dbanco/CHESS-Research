lambdaVals = logspace(-2,0,150);
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];


sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'dissertation_adjust2'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true','mcdl'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s0 = 1;
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

prefix = ['8_25_indep_',dataset];
topDir = ['E:\MCDLOF_processing\Outputs_',prefix,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

% criterion = 'discrepancy';
% criterion = 'truth_errorr';
% criterion = 'relaxed discrepancy';
% criterion = 'l-curve';
criterion = 'discrepancy range';
% criterion = 'triangle';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,3);
objectives = cell(NN,1);

useMin = 1;
relax_param = 1.15;
sig_ind = 2:6;
makeFigures = true;

for n = sig_ind
    inDir = [topDir,'\results_trial_1_sig_',num2str(n)];
    resultFile = [topDir,'\results_sig_',num2str(n),'.mat'];
    if exist(resultFile,'file')
        load(resultFile)
    else
        results = compute_metrics(inDir,sigmas(n),dataset,useMin);
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
    end
    objectives{n} = objective;
end
LcurveFile = fullfile(topDir,'l-curve_plot.png'); 
fig = figure(fig_num);
fig.Position = [100 100 1400 800];
saveas(gcf, LcurveFile);
removeWhiteSpace(LcurveFile)

%% Compute meanSNR
% [meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);
% 
%% Compute errors
% dirStartS = 'results_1';
% [~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(0,dataset);
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

%% Next copy figures associated with selected parameters to a folder
% pptFile = fullfile(topDir,['sim2_',prefix,...
%            '_relax_',num2str(relax_param),...
%            '_useMin0_',criterion,'_',opt.Penalty,'_',dirStartS,'.pptx']);
% titleStr = ['Sim 2 Param Select,',dirStartS];
% % createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)
% createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,sigmas,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind,objectives)