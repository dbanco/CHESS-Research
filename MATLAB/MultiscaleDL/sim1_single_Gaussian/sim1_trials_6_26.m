lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,19)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

fig_num = 22;

datasets = {'sim2_gaussian_tooth_matched','sim2_gaussian_tooth_unmatched',... %1,2
            'sim2_gaussian_tooth_matched2','sim2_gaussian_tooth_unmatched2',...%3,4
            'dissertation','dissertation_long',...%5,6
            'dissertation_long_separate','pseudo-voigt_unmatched',...%7,8
            'voigt_tooth_matched',...%9
            'gaussian_tooth_matched',...%10
            'gaussian_tooth_matched_long',...%11
            'gaussian_tooth_matched_long2'};%12
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

test_name = ['Outputs_6_26trials_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

test_name2 = ['Outputs_6_26trials_',dataset,'_',opt.Penalty,...
    '_D',opt.dictInit,num2str(opt.Dfixed),...
    '_X',opt.coefInit,num2str(opt.Xfixed),...
    '_recenter',num2str(opt.Recenter)];

topDir = ['E:\MCDLOF_processing\',test_name];
topDir2 = ['E:\MCDLOF_processing\',test_name2];

useMin = 1;
sig_ind = 1:6;
num_trials = 50;

objectives_indep = cell(numel(sig_ind),1);
objectives_of = cell(numel(sig_ind),1);
for n = sig_ind
    objective_indep = eval_trials(topDir,n,sigmas(n),dataset,useMin,num_trials,true);
    objectives_indep{n} = objective_indep;    
    objective_of = eval_trials(topDir2,n,sigmas(n),dataset,useMin,num_trials,false);
    objectives_of{n} = objective_of;
end

%% Plot Results
sig_ind = 1:6;
[meanSNR,noiseError] = computeSNR_noiseError(dataset,sig_ind);

error_stats_indep = compute_error_stats(objectives_indep,sig_ind);
error_stats_of = compute_error_stats(objectives_of,sig_ind);

% True Recon Error Figure
figure
hold on
errorbar(meanSNR,error_stats_indep.avg_true_error,...
                 error_stats_indep.std_true_error,'s-')
errorbar(meanSNR,error_stats_of.avg_true_error,...
                 error_stats_of.std_true_error,'s-')
title('Average over 50 trials')
xlabel('SNR','Fontsize',14)
ylabel('$\|\hat{{\bf b}}-{\bf f}\|_2$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')

% Dictionary Error Figure
figure
hold on
errorbar(meanSNR,error_stats_indep.avg_Derror,...
                 error_stats_indep.std_Derror,'s-')
errorbar(meanSNR,error_stats_of.avg_Derror,...
                 error_stats_of.std_Derror,'s-')
title('Average over 50 trials')
xlabel('SNR','Fontsize',14)
ylabel('$\|\hat{{\bf D}}-{\bf D}\|_2/\|{\bf D}\|_2$','Fontsize',14,...
       'interpreter','latex','Fontsize',14)
legend('MCDL',...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
       'MCDL-OF')
