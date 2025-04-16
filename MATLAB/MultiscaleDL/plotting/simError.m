function [noiseNorm,trueErrS,dataErrS,trueErrOF,dataErrOF] = simError(y_true,sigmas,sig_ind,topDir,dirStartS,selected_lam_s_vec,lambdaVals,dirStartOF,selected_lam_of_vec,lambdaOFVals)
%UNTITLED Summary of this function goes here

NN = numel(sigmas);

trueErrS = zeros(NN,1);
dataErrS = zeros(NN,1);
trueErrOF = zeros(NN,1);
dataErrOF = zeros(NN,1);
noiseNorm = zeros(NN,1);

if nargin < 8
    selected_lam_of_vec = zeros(NN,1);
    lambdaOFVals = [0,0,0];
    dirStartOF = 'none';
end

for n = sig_ind
    spDir = [dirStartS,'_sig_',num2str(n)];
    j_s = find(lambdaVals == selected_lam_s_vec(n));
    j_of = 1;
    j_hs = 1;
    dataFileS = sprintf("output_j%i_%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e_lam3_%0.2e.mat",...
                    j_s,j_of,j_hs,sigmas(n),selected_lam_s_vec(n),0,0);
    
    load(fullfile(topDir,spDir,dataFileS))
    D = outputs.D;
    N = outputs.N;
    M = outputs.M;
    y = outputs.y;
    X = outputs.X;
    scales = outputs.scales;
    center = (M+1)/2;
    
    AD = reSampleCustomArrayCenter(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);
    
    trueErrS(n) = norm(y_true-Yhat);
    dataErrS(n) = norm(squeeze(y)-Yhat);
    
    if nargin > 7
        hsDir = [dirStartOF,'_sig_',num2str(n)];
    
        j_of = find(lambdaOFVals == selected_lam_of_vec(n));
        dataFileOF = sprintf("output_j%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e",...
                            j_of,sigmas(n),selected_lam_s_vec(n),selected_lam_of_vec(n));

        load(fullfile(topDir,hsDir,dataFileOF))
        D = outputs.D;
        N = outputs.N;
        M = outputs.M;
        y = outputs.y;
        X = outputs.X;
        scales = outputs.scales;
        center = (M+1)/2;
        
        AD = reSampleCustomArrayCenter(N,D,scales,center);
        AD = padarray(AD,[0 M-1 0],0,'post');
        ADf = fft2(AD);
        Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
        Yhat = gather(Yhat);
        
        trueErrOF(n) = norm(y_true-Yhat);
        dataErrOF(n) = norm(squeeze(y)-Yhat);
    end
    
    T = outputs.T;
    noiseNorm(n) = norm(randn(N,T)*sigmas(n));
%     noiseNorm2(n) = norm(y_true-squeeze(y));

    if 0
        ff = figure();
        hold on
        kk = 1;
        for ttt = [1,11,21]
            subplot(3,1,kk)
            plot(outputs.y(:,:,ttt),'-')
            hold on
            plot(y_true(:,ttt),'-')
            plot(Yhat(:,ttt),'-')
            kk = kk + 1; 
        end
        ff.Position = [941 59 560 825];
        legend('data','truth','recon')
    end

end

trueErrS = trueErrS(sig_ind);
dataErrS = dataErrS(sig_ind);
trueErrOF = trueErrOF(sig_ind);
dataErrOF = dataErrOF(sig_ind);
noiseNorm = noiseNorm(sig_ind);

end