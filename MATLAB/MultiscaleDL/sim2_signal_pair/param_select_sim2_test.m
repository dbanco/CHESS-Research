% lambdaVals = logspace(-3,1,120);
% lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,5)];
% lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1];

lambdaVals = logspace(-3,1,120);
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,5)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1];

sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'sim2_gaussian_tooth_matched','sim2_gaussian_tooth_unmatched',...
            'sim2_gaussian_tooth_matched2','sim2_gaussian_tooth_unmatched2',...
            'dissertation','dissertation_long','dissertation_long_separate',...
            'pseudo-voigt_unmatched','voigt_tooth_matched',...
            'gaussian_tooth_matched','gaussian_tooth_matched_long',...
            'gaussian_tooth_matched_long2','dissertation_adjust2'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s0 = 13;
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

test_name = ['Outputs_7_16_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\',test_name];

topDir2 = ['E:\MCDLOF_processing\',test_name];

figDir = [topDir,'_sig_',num2str(i)];

% criterion = 'discrepancy';
% criterion = 'truth_error';
% criterion = 'relaxed discrepancy';
% criterion = 'discrepancy range';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,3);
objectives = cell(NN,1);

useMin = 1;
relax_param = 1.5;
sig_ind = 2:6;
makeFigures = 1;
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
        psFigDir = fullfile(topDir,['ps_',criterion,'_',num2str(relax_param)]);
        generateFiguresToy1zpad_center(psFigDir,outputs,suffix,[4,8],useMin);
        close all
    end
    objectives{n} = objective;
end

%% Select minimizers
% nn = find(min(results.true_error) == results.true_error);
% lam_vals = results.lambda_vec(nn,:)
% lam_inds = convertToInds(lam_vals,lambdaVals,lambdaOFVals,lambdaHSVals)
% 
% nn = find(min(results.X_error) == results.X_error);
% lam_vals = results.lambda_vec(nn,:)
% lam_inds = convertToInds(lam_vals,lambdaVals,lambdaOFVals,lambdaHSVals)
% 
% nn = find(min(results.vdf_error) == results.vdf_error);
% lam_vals = results.lambda_vec(nn,:)
% lam_inds = convertToInds(lam_vals,lambdaVals,lambdaOFVals,lambdaHSVals)
% 
% nn = find(min(results.shift_error) == results.shift_error);
% lam_vals = results.lambda_vec(nn,:)
% lam_inds = convertToInds(lam_vals,lambdaVals,lambdaOFVals,lambdaHSVals)
% 
% nn = find(min(results.l0_norm) == results.l0_norm);
% lam_vals = results.lambda_vec(nn,:)
% lam_inds = convertToInds(lam_vals,lambdaVals,lambdaOFVals,lambdaHSVals)


results.lambda_vec(65,:)


%% Compute meanSNR
% [meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);

%% Compute errors
% dirStartS = 'results_trial_1';
% % selected_lam_all_vec2 = [43,1,1;44,1,1;55,1,1;58,1,1;65,1,1;69,1,1];
% [~,y_true,~,~,~] = sim_switch_multiscale_dl(sigmas(n),dataset);
% [noiseNorm,trueErr1,dataErr1,l0_norm1,trueErr2,dataErr2,l0_norm2] = simError(y_true,sigmas,sig_ind,topDir,dirStartS,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals);
% 
% figure()
% hold on
% plot(meanSNR,noiseNorm,'s-')
% plot(meanSNR,trueErr2,'o-')
% plot(meanSNR,trueErr1,'o-')
% xlabel('SNR','Fontsize',14)
% ylabel('Error','Fontsize',14)
% legend('$\|{\bf w}\|_2$','$\|\hat{{\bf b}}-{\bf f}\|_2$',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%     '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...
%     'interpreter','latex','Fontsize',14)
% 
% figure()
% hold on
% plot(meanSNR,l0_norm2,'o-')
% plot(meanSNR,l0_norm1,'o-')
% xlabel('SNR','Fontsize',14)
% ylabel('$l_0$-norm','Fontsize',14)
% legend('$\|{\bf{X}\|_0$',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%     '$\|{\bf{X}\|_0$ (OF)',...
%     'interpreter','latex','Fontsize',14)
% 
% 
% %% Next copy figures associated with selected parameters to a folder
% pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim2_of_',test_name,'_relax_',num2str(relax_param),'.pptx'];
% titleStr = ['Sim 2 OF/HS Param Select,',dirStartS];
% % createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)
% createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,sigmas,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind,objectives)