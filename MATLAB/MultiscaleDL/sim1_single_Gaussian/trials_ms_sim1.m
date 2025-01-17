lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-4 5e-4 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];

sigmas = 0:0.01:0.1;
NN = numel(sigmas)-1;

dataset = 'steps_matched';
penalty = 'log';
testType = 'Dflat0_Xzeros0';
topDir = 'E:\Outputs_sim1_trials\';

num_trials = 19;
[true_error_s,...
 ~,~,~,~,~] = sim1_trials_results(sigmas,dataset,...
                                                    penalty,testType,...
                                                    topDir,num_trials);
[meanSNR,noiseError] = computeSNR_noiseError(dataset);

%% True and Data Error Averaged
figure(1)
hold on
x = meanSNR(2:9);
y1 = noiseError(2:9);
mu = mean(true_error(2:9,:),2);
sig = std(true_error(2:9,:)');

blue = [0 0.4470 0.7410];
errorbar(x,mu,sig,'Color',blue)
plot(x,mu,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',4)
plot(x,y1,'-o','MarkerSize',4)

ylabel('$\frac{1}{2}||{\bf f} - \hat{\bf f} ||_2$','interpreter','latex',...
    'FontSize',16)
xlabel('SNR','FontSize',16)
legend('Noise Error','MCDL')

