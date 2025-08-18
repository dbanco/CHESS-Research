% lambdaVals = logspace(-3,1,120);
% lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,5)];
% lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1];

lambdaVals = logspace(-3,1,120);
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,5)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1];

sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'mmpad_ring3'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s_data = 1;
s_pen = 2;
s_xinit = 1;
s_dinit = 2;
s_dfix = 1;
s_recenter = 1;
dataset = datasets{s_data};
opt.Penalty = penalties{s_pen};
opt.coefInit = xinits{s_xinit};
opt.dictInit = dinits{s_dinit};
opt.Dfixed = dfixes{s_dfix};
opt.Recenter = recenters{s_recenter};
opt.Xfixed = 0;


for K = 1:4
for M = [51,101,151,201]

opt.M = M;
n = 1;

test_name = ['Outputs_8_14_',dataset,'_',opt.Penalty,...
    '_M',num2str(opt.M),...
    '_K',num2str(K),...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\',test_name];
topDir2 = ['E:\MCDLOF_processing\',test_name];
figDir = [topDir,'_M_',num2str(opt.M),'_K_',num2str(K)];

% criterion = 'discrepancy';
% criterion = 'truth_error';
% criterion = 'relaxed discrepancy';
criterion = 'discrepancy range';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,3);
objectives = cell(NN,1);

useMin = 0;
relax_param = 1.1;
makeFigures = 1;

[y,~,~,~,~,~,~] = sim_switch_multiscale_dl(0,dataset);
noise_est = sqrt(estimate_noise(y,5));

inDir = [topDir,'\results_K_',num2str(K)];
resultFile = [topDir,'\results_K_',num2str(K),'.mat'];
if exist(resultFile,'file')
    load(resultFile)
else
    results = compute_metrics(inDir,0,dataset,useMin,true);
    save(resultFile,'results')
end
[lambda_all,objective] = param_select_mcdl(results,criterion,noise_est,dataset,relax_param,fig_num);
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
                  j_s,j_of,j_hs,0,outputs.lambda,outputs.lambda2,outputs.lambda3);
    psFigDir = fullfile(topDir,['ps_',criterion,'_',num2str(relax_param),...
                                      '_useMin',num2str(useMin)]);
    generateFiguresToy1zpad_center(psFigDir,outputs,suffix,[4,8],useMin);
    close all
end
end
end
%objectives{n} = objective;
