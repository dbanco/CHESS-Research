lambdaVals = logspace(-3,1,120);
% lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,5)];
% lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1];
lambdaOFVals = [0 logspace(-4,0,20)];
lambdaHSVals = [0 logspace(-4,0,10)];

sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'dissertation_adjust2'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
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

test_name = ['Outputs_7_23a1_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

test_name2 = ['Outputs_8_18of_',dataset,'_',opt.Penalty,...
    '_D','mcdl',num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\',test_name];
topDir = ['E:\MCDLOF_processing\',test_name2];

% criterion = 'discrepancy';
% criterion = 'truth_errorr';
% criterion = 'relaxed discrepancy';
criterion = 'discrepancy range';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,3);
objectives = cell(NN,1);

useMin = 0;
relax_param = 1.1;
sig_ind = 2:6;
makeFigures = true;

for n = sig_ind
    inDir = [topDir,'\results_trial_1_sig_',num2str(n)];
    [lambda_all,objective] = param_select_3D(inDir,fig_num,criterion,sigmas(n),dataset,useMin,relax_param,false);
    selected_lam_all_vec(n,:) = lambda_all;
    selected_lam_s_vec(n) = lambda_all(1);
    selected_inds(n,1) = find(selected_lam_s_vec(n) == lambdaVals);
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
fig = gcf;
fig.Position = [100 100 1400 800];
saveas(gcf, LcurveFile);
removeWhiteSpace(LcurveFile)

%% Compute meanSNR
[meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);

%% Compute errors
dirStartS = 'results_trial_1';
% selected_lam_all_vec2 = [0,0,0;9,1,1;12,1,1;21,1,1;30,1,1;37,1,1];
[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(0,dataset);
[noiseNorm,trueErr1,dataErr1,l0_norm1,trueErr2,dataErr2,l0_norm2] = simError(y_true,sigmas,sig_ind,topDir,dirStartS,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals);

figure()
hold on
plot(meanSNR,noiseNorm,'s-')
plot(meanSNR,trueErr2,'o-')
plot(meanSNR,trueErr1,'o-')
xlabel('SNR','Fontsize',14)
ylabel('Error','Fontsize',14)
legend('$\|{\bf w}\|_2$','$\|\hat{{\bf b}}-{\bf f}\|_2$',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...
    'interpreter','latex','Fontsize',14)

% figure()
% hold on
% plot(meanSNR,l0_norm2,'o-')
% plot(meanSNR,l0_norm1,'o-')
% xlabel('SNR','Fontsize',14)
% ylabel('$l_0$-norm','Fontsize',14)
% legend('$\|{\bf{X}\|_0$',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%     '$\|{\bf{X}\|_0$ (OF)',...
%     'interpreter','latex','Fontsize',14)


%% Next copy figures associated with selected parameters to a folder
% pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim1_of_',test_name,'_relax_',num2str(relax_param),'.pptx'];
% titleStr = ['Sim 1 Param Select,',dirStartS];
% createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)
% createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,sigmas,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind,objectives)