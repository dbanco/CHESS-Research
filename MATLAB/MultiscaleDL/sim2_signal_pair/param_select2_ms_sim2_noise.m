lambdaVals = [1e-4,5e-4,1e-3,5e-3,1e-2,2e-2,linspace(3e-2,8e-1,100)];
lambdaHSVals = [0 1e-3 2e-3];
lambdaOFVals = [0 1e-4,5e-4,1e-3,5e-3,1e-2,linspace(5e-2,1,50)];
sigmas = 0:0.01:0.1;
NN = numel(sigmas);

tradeoff_s = 0.75;
tradeoff_of = 1;
scaleP = [0.157,5.12,19.02,52.406,0,100];

for j_hs = 3 %2-3
    selected_lam_s_vec = zeros(NN,1);
    selected_lam_of_vec = zeros(NN,1);
    
    % close all
    topDir = 'C:\Users\dpqb1\Documents\Outputs2024_8_29\';
    fig_num = 12;
    for n = 2:NN
        inDir = [topDir,'\signal_pair_8_29_24_X0_D0_V00_sig_',num2str(n)];
        [lambda_s_sel,j_s] = param_select_lambda_s(inDir,tradeoff_s,scaleP,11);
        selected_lam_s_vec(n) = lambda_s_sel;
    
        outDir = ['signal_pair_8_29_24_coupled',num2str(lambdaHSVals(j_hs)),'_sig_',num2str(n)];
        folderPath = fullfile(topDir,outDir);
        [lambda_of_sel,j_of] = param_select_lambda_of(folderPath,tradeoff_of,scaleP,fig_num);
        selected_lam_of_vec(n) = lambda_of_sel;
    end
    
    %% Next plot the true recovery as a function of SNR
    sigmas = 0:0.01:0.1;
    meanSNR = zeros(numel(sigmas),1);
    for n = 2:numel(sigmas)
        [y,y_true,N,M,T] = generate_signal_pair(sigmas(n));
        
        SNR = zeros(T,1);
        nPwr = zeros(T,1);
        sigPwr = zeros(T,1);
        for t = 1:T
            SNR(t) = norm(y_true(:,t))/norm(y(:,t)-y_true(:,t));
            nPwr(t) = norm(y(:,t)-y_true(:,t));
            sigPwr(t) = norm(y_true(:,t));
        end
        
        meanSNR(n) = mean(SNR);
    end
    %% Need to write a function to compute the various errors for all noise levels
    % and that maybe just takes the file name as input so that it works on both
    
    outDirOF = "C:\Users\dpqb1\Documents\Outputs2024_8_29\";
    dirStartS = 'signal_pair_8_29_24_X0_D0_V00';
    dirStartOF = ['signal_pair_8_29_24_coupled',num2str(lambdaHSVals(j_hs))];
       
    dataFile = sprintf("output_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e",...
                        j_of,sigmas(n),selected_lam_s_vec(n),selected_lam_of_vec(n));
    
    [noiseNorm,trueErrS,dataErrS,trueErrOF,dataErrOF] = simError(sigmas,outDirOF,dirStartS,dirStartOF,...
                                selected_lam_s_vec,selected_lam_of_vec,lambdaOFVals,y_true);
        
    %%
    figure()
    hold on
    plot(meanSNR(2:NN),noiseNorm(2:NN),'s-')
    plot(meanSNR(2:NN),trueErrS(2:NN),'o-')
    % plot(meanSNR(2:NN),dataErrS(2:NN),'x-')
    plot(meanSNR(2:NN),trueErrOF(2:NN),'o-')
    % plot(meanSNR(2:NN),dataErrOF(2:NN),'x-')
    % plot(meanSNR(2:NN),noiseNorm2(2:NN),'x-')
    xlabel('SNR','Fontsize',14)
    ylabel('Error','Fontsize',14)
    legend('$\|{\bf w}\|_2$','$\|\hat{{\bf b}}-{\bf f}\|_2$',...
        '$\|\hat{{\bf b}}-{\bf f}\|_2$ (OF)',...
        'interpreter','latex','Fontsize',14)
end
%% Next copy figures associated with selected parameters to a folder

% pptFile = 'C:\Users\dpqb1\Documents\MCDL Paper\recons_dicts_sim2.pptx';
% titleStr = 'Sim 2 Recovery';
% createPowerpointSim(pptFile,titleStr,topDir,sigmas,dirStartS,dirStartOF,...
%     lambdaOFVals,lambdaHSVals,selected_lam_of_vec,selected_lam_s_vec,meanSNR)


pptFile = 'C:\Users\dpqb1\Documents\MCDL Paper\recons_dicts_sim2_indep_0.75.pptx';
titleStr = 'Sim 2 Recovery';
createPowerpointSimS(pptFile,titleStr,topDir,sigmas,dirStartS,dirStartOF,...
    selected_lam_s_vec,meanSNR)

