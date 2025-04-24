lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'sim2_gaussian_tooth_matched','sim2_gaussian_tooth_unmatched',...
            'sim2_gaussian_tooth_matched2','sim2_gaussian_tooth_unmatched2',...
            'dissertation','dissertation_long','dissertation_long_separate',...
            'pseudo-voigt_unmatched','voigt_tooth_matched'};
penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};
recenters = {0,1};

% Setup Dataset
s0 = 9;
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
topDir = ['E:\MCDLOF_processing\Outputs_4_23',...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

% criterion = 'discrepancy';
% criterion = 'truth_error';
criterion = 'relaxed discrepancy';

selected_lam_s_vec = zeros(NN,1);
selected_lam_of_vec = zeros(NN,1);
selected_inds = zeros(NN,1);
sig_ind = 2:6;
for n = sig_ind
    inDir = [topDir,'\results_sig_',num2str(n)];
    [~,y_true,~,~,~] = sim_switch_multiscale_dl(sigmas(n),dataset);
    [lambda_s_sel,j_s] = param_select_3D(inDir,fig_num,criterion,sigmas(n),y_true,useMin);
    selected_lam_s_vec(n) = lambda_all(1);
    selected_inds(n) = find(selected_lam_s_vec(n) == lambdaVals);
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
pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim2_lam_s',criterion,'_',dataset,'_',dirStartS,'.pptx'];
titleStr = ['Sim 2 Recovery,',dirStartS];
createPowerpointSimS(pptFile,titleStr,meanSNR,topDir,sigmas,dirStartS,selected_lam_s_vec,lambdaVals,LcurveFile,criterion,sig_ind)

