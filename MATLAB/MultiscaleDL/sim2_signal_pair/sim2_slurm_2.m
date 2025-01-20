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
% opt.Penalty = 'l1-norm';
opt.Penalty = 'log';
opt.coefInit = 'zeros';
opt.dictInit = 'rand';
opt.a = 1;

% Multiscale dictionary setup
K = 2;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,8, 1/8);
scales{2} = genRationals([0;1],[1;1],16,8, 1/8);

penalties = {'l1-norm','log'};
xinits = {'zeros','true'};
dinits = {'rand','flat','true'};
dfixes = {0,1};

scriptFileName = 'mcdlof_bash.sh';
funcName = 'sim_mcdlof_wrapper';
jobDir = '/cluster/home/dbanco02/jobs/';
k = 1;

for r_ind = 1
for s1 = 2
    for s2 = 1
        for s3 = 1
            for s4 = 1

opt.Penalty = penalties{s1};
opt.coefInit = xinits{s2};
opt.dictInit = dinits{s3};
opt.Dfixed = dfixes{s4};

if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
    continue
end
for sig_i = 2:4
    % j_s_select = find(lambdaVals == selected_lam_s_vec(sig_i));
    for j_s = 1:10
        for j_hs = [2,3,4,5]  
            for j_of = 1:10:56
                dataset = 'sim2_gaussian_tooth_matched';
                topDir = ['/cluster/home/dbanco02/Outputs_sim2_lam_all',num2str(j_hs),...
                    '_D',opt.dictInit,num2str(opt.Dfixed),...
                    '_X',opt.coefInit,num2str(opt.Xfixed),'/',...
                    dataset,'_',opt.Penalty,'_results'];
                
                varin = {lambdaVals,lambdaOFVals,lambdaHSVals,...
                        j_s,j_of,j_hs,sigmas,sig_i,opt,topDir,dataset,K,scales};
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