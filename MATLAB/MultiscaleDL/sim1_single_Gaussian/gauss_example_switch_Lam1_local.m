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
opt.coefInit = 'zeros';
opt.dictInit = 'flat';
opt.maxLayers = 2;
opt.delta = 0.02;
opt.epsilon = 0.05;

k = 1;
scriptFileName = 'mcdlof_bash.sh';
jobDir = '/cluster/home/dbanco02/jobs/';
for i = [1]
    for j_s = [1,5,10,100]
        for j_hs = 1  
            for j_of = 1
                % topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_23_Dflat_Xzero\unmatched_results';
                % figDir = [topDir,'_sig_',num2str(i)];
                % mkdir(figDir)
                % dataset = 'unmatched';
                % mcdlof_wrapper_sim1_switch(lambdaVals,lambdaOFVals,lambdaHSVals,...
                %                             j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir,dataset);
                % 
                % topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_23_Dflat_Xzero\matched_results';
                % figDir = [topDir,'_sig_',num2str(i)];
                % mkdir(figDir)
                % dataset = 'matched';
                % mcdlof_wrapper_sim1_switch(lambdaVals,lambdaOFVals,lambdaHSVals,...
                %                             j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir,dataset);

                % topDir = 'C:\Users\dpqb1\Documents\Outputs2024_10_24_Dtrue_Xzero\steps_results';
                % figDir = [topDir,'_sig_',num2str(i)];
                % mkdir(figDir)
                % dataset = 'steps';
                % mcdlof_wrapper_sim1_switch(lambdaVals,lambdaOFVals,lambdaHSVals,...
                %                             j_s,j_of,j_hs,sigmas,i,opt,K,scales,topDir,dataset);
                dataset = 'steps_matched';
                topDir = ['C:\Users\dpqb1\Documents\Outputs2024_11_14',...
                    '_D',opt.dictInit,num2str(opt.Dfixed),...
                    '_X',opt.coefInit,num2str(opt.Xfixed),'\',dataset,'_results'];
                 
                mcdlof_wrapper_sim1_switch(lambdaVals,lambdaOFVals,lambdaHSVals,...
                                            j_s,j_of,j_hs,sigmas,i,opt,topDir,dataset);
            end
        end
    end
end

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