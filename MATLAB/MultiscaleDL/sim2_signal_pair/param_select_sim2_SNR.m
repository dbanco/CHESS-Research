lambdaVals = logspace(-2,0,150);
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];

dataset = 'gaus_loren_nooverlap';
SNRs = [20,16,12,8,4];
sigmas = zeros(numel(SNRs),1);
sig_ind = 1:5;
for i = sig_ind
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
end

NN = numel(sigmas);

fig_num = 22;

penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true','mcdl'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s = [1,2,1,2,1,1];
opt.Penalty = penalties{s(2)};
opt.coefInit = xinits{s(3)};
opt.dictInit = dinits{s(4)};
opt.Dfixed = dfixes{s(5)};
opt.Recenter = recenters{s(6)};
opt.Xfixed = 0;

prefix = ['9_8_of_',dataset];
topDir = ['E:\MCDLOF_processing\Outputs_',prefix,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

% criterion = 'discrepancy';
% criterion = 'truth_errorr';
% criterion = 'relaxed discrepancy';
% criterion = 'l-curve';
criterion = 'discrepancy range';
% criterion = 'discrepancy range of-log';
% criterion = 'triangle';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,3);
objectives = cell(NN,1);

useMin = 1;
relax_param = 1.1;
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