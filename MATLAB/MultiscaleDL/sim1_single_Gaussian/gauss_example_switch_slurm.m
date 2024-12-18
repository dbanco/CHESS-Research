%% Multiscale 1D dictionary learning toy problem
% Directory
lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];

% Experiment Setup
sigmas = 0:0.01:0.1;

% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 1000;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho and sigma params
opt.rho = 300;
opt.sigma = 50;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.UpdateVelocity = 1;
opt.HSiters = 100;
opt.useGpu = 0;
opt.Xfixed = 0;
opt.Dfixed = 0;
opt.Penalty = 'l1-norm';
% opt.Penalty = 'log';
opt.coefInit = 'zeros';
opt.dictInit = 'flat';
opt.a = 1;

penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'flat','true'};
dfixes = {0,1};

load('/cluster/home/dbanco02/CHESS-Research/MATLAB/MultiscaleDL/sim1_single_Gaussian/sim1_selected_lam_s.mat')

scriptFileName = 'mcdlof_bash.sh';
funcName = 'mcdlof_wrapper_sim1_switch';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

for r_ind = 1
for s1 = 1
    for s2 = 2
        for s3 = 1
            for s4 = 1

opt.Dfixed = dfixes{s1};
opt.Penalty = penalties{s2};
opt.coefInit = xinits{s3};
opt.dictInit = dinits{s4};
if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
    continue
end
for sig_i = 3
    j_s_select = find(lambdaVals == selected_lam_s_vec(sig_i));
    for j_s = j_s_select
        for j_hs = 2:4  
            for j_of = 1:56
                dataset = 'steps_matched';
                topDir = ['/cluster/home/dbanco02/Outputs_lam_of',num2str(r_ind),...
                    '_D',opt.dictInit,num2str(opt.Dfixed),...
                    '_X',opt.coefInit,num2str(opt.Xfixed),'/',...
                    dataset,'_',opt.Penalty,'_results'];
                
                varin = {lambdaVals,lambdaOFVals,lambdaHSVals,...
                        j_s,j_of,j_hs,sigmas,sig_i,opt,topDir,dataset};
                save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
                k = k + 1;
            end
        end
    end
end
            end
        end
    end
end
end

slurm_write_bash(k-1,jobDir,scriptFileName,sprintf('1-%i',k-1))
%% Compare data
%
% [yn1,y_m,N,M,T,~,~] = gaus_example_matched_multiscale_dl(0);
% [yn2,y_um,N,M,T,~,~] = gaus_example_unmatched_multiscale_dl(0);
% figure(2)
% subplot(3,1,1)
% imagesc(y_m)
% clim([0,0.4])
% colorbar()
% 
% subplot(3,1,2)
% imagesc(y_um)
% clim([0,0.4])
% colorbar()
% 
% subplot(3,1,3)
% imagesc(abs(y_m -y_um))
% colorbar()
% 
% figure(3)
% plot(y_um(:,10))
% hold on
% plot(y_m(:,10))