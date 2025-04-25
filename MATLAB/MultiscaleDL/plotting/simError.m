function [noiseNorm,trueErrS,dataErrS,l0_normS,trueErrOF,dataErrOF,l0_normOF] = simError(y_true,sigmas,sig_ind,topDir,dirStart1,selected_lam_all_vec,lambdaVals,lambdaOFVals,lambdaHSVals,topDir2,selected_lam_all_ind_vec2)
%UNTITLED Summary of this function goes here

NN = numel(sigmas);

trueErrS = zeros(NN,1);
dataErrS = zeros(NN,1);
trueErrOF = zeros(NN,1);
dataErrOF = zeros(NN,1);
noiseNorm = zeros(NN,1);
l0_normS = zeros(NN,1);
l0_normOF = zeros(NN,1);

if nargin < 10
    selected_lam_all_vec2 = zeros(NN,1);
end

for n = sig_ind
    spDir = [dirStart1,'_sig_',num2str(n)];
    j_s = find(lambdaVals == selected_lam_all_vec(n,1));
    j_of = find(lambdaOFVals == selected_lam_all_vec(n,2));
    j_hs = find(lambdaHSVals == selected_lam_all_vec(n,3));
    dataFileS = sprintf("output_j%i_%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e_lam3_%0.2e.mat",...
                    j_s,j_of,j_hs,sigmas(n),...
                    selected_lam_all_vec(n,1),...
                    selected_lam_all_vec(n,2),...
                    selected_lam_all_vec(n,3));
    
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
    l0_normS(n) = sum(X(:)>0,'all')
    
    if nargin > 10
        hsDir = [dirStart1,'_sig_',num2str(n)];

        j_s = selected_lam_all_ind_vec2(n,1)
        j_of = selected_lam_all_ind_vec2(n,2)
        j_hs = selected_lam_all_ind_vec2(n,3)
        dataFileOF = sprintf("output_j%i_%i_%i_sig_%0.2e_lam1_%0.2e_lam2_%0.2e_lam3_%0.2e.mat",...
                        j_s,j_of,j_hs,sigmas(n),...
                        lambdaVals(j_s),...
                        lambdaOFVals(j_of),...
                        lambdaHSVals(j_hs));
    
        load(fullfile(topDir2,hsDir,dataFileOF))
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
        l0_normOF(n) = sum(X(:)>0,'all')
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
l0_normS = l0_normS(sig_ind);
trueErrOF = trueErrOF(sig_ind);
dataErrOF = dataErrOF(sig_ind);
l0_normOF = l0_normOF(sig_ind);
noiseNorm = noiseNorm(sig_ind);

end