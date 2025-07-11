function [SNR,noiseError,noiseTheor] = computeSNR_noiseError(dataset,sig_ind)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
sigmas = 0:0.01:0.1;
NN = numel(sigmas);
sigPwr = zeros(numel(sigmas),1);
noisePwr = zeros(numel(sigmas),1);
SNR = zeros(numel(sigmas),1);
noiseError = zeros(numel(sigmas),1);
noiseTheor = zeros(numel(sigmas),1);
for n = sig_ind
    [y,y_true,N,M,T] = sim_switch_multiscale_dl(sigmas(n),dataset);
    
    % SNR_time = zeros(T,1);
    % nPwr = zeros(T,1);
    % sigPwr = zeros(T,1);
    % for t = 1:T
    %     SNR_time(t) = norm(y_true(:,t))^2/norm(y(:,t)-y_true(:,t))^2;
    %     nPwr(t) = norm(y(:,t)-y_true(:,t))^2;
    %     sigPwr(t) = norm(y_true(:,t))^2;
    % end
    sigPwr(n) = mean(y(:).^2);
    noisePwr(n) = mean((y(:)-y_true(:)).^2);
    SNR(n) = 10*log10(sigPwr(n)/noisePwr(n));
    noiseTheor(n) = sigmas(n)^2;
end

SNR = SNR(sig_ind);
noiseError = noiseError(sig_ind);
noiseTheor = noiseTheor(sig_ind);

end