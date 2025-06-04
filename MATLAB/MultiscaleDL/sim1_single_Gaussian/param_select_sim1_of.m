lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
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
s5 = 2;
dataset = datasets{s0};
opt.Penalty = penalties{s1};
opt.coefInit = xinits{s2};
opt.dictInit = dinits{s3};
opt.Dfixed = dfixes{s4};
opt.Recenter = recenters{s5};
opt.Xfixed = 0;

topDir = ['E:\MCDLOF_processing\Outputs_4_24_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir2 = ['E:\MCDLOF_processing\Outputs_4_19_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

% criterion = 'discrepancy';
% criterion = 'truth_errorr';
% criterion = 'relaxed discrepancy';
criterion = 'discrepancy range';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_lam_all_vec = zeros(NN,3);
selected_inds = zeros(NN,1);
objectives = cell(NN,1);

useMin = 1;
relax_param = 1.05;
sig_ind = 2:4;
for n = sig_ind
    inDir = [topDir,'\results_sig_',num2str(n)];
    [lambda_all,objective] = param_select_3D(inDir,fig_num,criterion,sigmas(n),dataset,useMin,relax_param);
    selected_lam_s_vec(n) = lambda_all(1);
    selected_lam_all_vec(n,:) = lambda_all;
    selected_inds(n) = find(selected_lam_s_vec(n) == lambdaVals);
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
dirStartS = 'results';
selected_lam_all_vec2 = [0,0,0;9,1,1;12,1,1;21,1,1;30,1,1;37,1,1];
[~,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(0,dataset);
[noiseNorm,trueErr1,dataErr1,l0_norm1,trueErr2,dataErr2,l0_norm2] = simError(y_true,sigmas,sig_ind,topDir,dirStartS,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,topDir2,selected_lam_all_vec2);

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

figure()
hold on
plot(meanSNR,l0_norm2,'o-')
plot(meanSNR,l0_norm1,'o-')
xlabel('SNR','Fontsize',14)
ylabel('$l_0$-norm','Fontsize',14)
legend('$\|{\bf{X}\|_0$',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    '$\|{\bf{X}\|_0$ (OF)',...
    'interpreter','latex','Fontsize',14)


%% Next copy figures associated with selected parameters to a folder
pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim1_lam_of_4_24',...
           '_relax_',num2str(relax_param),...
           '_useMin1_',criterion,'_',opt.Penalty,'_',dirStartS,'.pptx'];
titleStr = ['Sim 1 Param Select,',dirStartS];
% createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)
createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,sigmas,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind,objectives)