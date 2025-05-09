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

count = 0;
total = 0;

for s1 = 1
    for s2 = 1:2
        for s3 = 1
            for s4 = 1


opt.Dfixed = dfixes{s1};
opt.Penalty = penalties{s2};
opt.coefInit = xinits{s3};
opt.dictInit = dinits{s4};
if (opt.Dfixed == 1) && strcmp(opt.dictInit, 'flat')
    continue
end
k = 1;
for sig_i = 2
    for j_s = [1:10,16:5:106]
        for j_hs = 1  
            for j_of = 1
                dataset = 'steps_matched';
                topDir = ['E:\Outputs2024_12_3',...
                    '_D',opt.dictInit,num2str(opt.Dfixed),...
                    '_X',opt.coefInit,num2str(opt.Xfixed),'\',...
                    dataset,'_',opt.Penalty,'_results'];
                
                imageFile = dir(fullfile([topDir,'_sig_',num2str(sig_i)],...
                                         ['vdf_j', num2str(j_s), '_*.png']));
                total = total + 1
                if numel(imageFile) == 1
                    continue
                end
                count = count + 1
                mcdlof_wrapper_sim1_switch(lambdaVals,lambdaOFVals,lambdaHSVals,...
                                            j_s,j_of,j_hs,sigmas,sig_i,opt,topDir,dataset);
            end
        end
    end
end
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