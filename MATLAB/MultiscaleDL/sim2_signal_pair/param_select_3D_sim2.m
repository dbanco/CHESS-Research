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

suffix = [dataset,'_',opt.Penalty,...
        '_D',opt.dictInit,num2str(opt.Dfixed),...
        '_X',opt.coefInit,num2str(opt.Xfixed),...
        '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\Outputs_4_18_',suffix];
% topDir = ['E:\MCDLOF_processing\Outputs_a0.1_',suffix];

criterion = 'discrepancy';
% criterion = 'truth_error';
useMin = true;

selected_lam_all_vec = zeros(NN,3);
selected_inds_all = zeros(NN,3);

sig_ind = 2:6;
true_error_sig = zeros(numel(sig_ind),1);
i = 1;
for n = sig_ind
    inDir = [topDir,'\results','_sig_',num2str(n)];
    [~,y_true,~,~,~] = sim_switch_multiscale_dl(sigmas(n),dataset);
    [lambda_all,objective] = param_select_3D(inDir,fig_num,criterion,sigmas(n),y_true,useMin);
    selected_lam_all_vec(n,:) = lambda_all;
    selected_inds_all(n,:) = lambdasToInds(lambda_all,{lambdaVals,lambdaOFVals,lambdaHSVals});
    true_error_sig(i) = objective.true_error;
    i = i + 1;
end
LcurveFile = fullfile(topDir,'l-curve_plot.png'); 
fig = gcf;
fig.Position = [100 100 1400 800];
saveas(gcf, LcurveFile);
removeWhiteSpace(LcurveFile)

[meanSNR,noiseError,noiseTheor] = computeSNR_noiseError(dataset,sig_ind);

figure()
hold on
plot(meanSNR,noiseTheor,'o-')
plot(meanSNR,true_error_sig,'s-')
xlabel('SNR','Fontsize',14)
ylabel('Error','Fontsize',14)
legend('$\sigma\sqrt{NT}$','$\|\hat{{\bf b}}-{\bf f}\|_2$',... 
        'interpreter','latex','Fontsize',14)
    % '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...



%% Next copy figures associated with selected parameters to a folder
pptFile = ['C:\Users\dpqb1\Documents\MCDL Paper\sim2_a1_lam_s_',criterion,'_',suffix,'.pptx'];
titleStr = ['Sim 2 Recovery, ',suffix];
createPowerpointSimAll(pptFile,titleStr,meanSNR,topDir,sigmas,...
    selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,LcurveFile,criterion,sig_ind)

