function [meanSNR,noiseError,noiseTheor] = computeSNR_noiseError(dataset,sig_ind)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
sigmas = 0:0.01:0.1;
NN = numel(sigmas);
meanSNR = zeros(numel(sigmas),1);
noiseError = zeros(numel(sigmas),1);
noiseTheor = zeros(numel(sigmas),1);
for n = sig_ind
    [y,y_true,N,M,T] = sim_switch_multiscale_dl(sigmas(n),dataset);
    
    SNR = zeros(T,1);
    nPwr = zeros(T,1);
    sigPwr = zeros(T,1);
    for t = 1:T
        SNR(t) = norm(y_true(:,t))/norm(y(:,t)-y_true(:,t));
        nPwr(t) = norm(y(:,t)-y_true(:,t));
        sigPwr(t) = norm(y_true(:,t));
    end
    
    noiseError(n) = norm(y(:)-y_true(:)); 
    meanSNR(n) = mean(SNR);
    noiseTheor(n) = sqrt(sigmas(n)^2*N*T);
end

meanSNR = meanSNR(sig_ind);
noiseError = noiseError(sig_ind);
noiseTheor = noiseTheor(sig_ind);

end