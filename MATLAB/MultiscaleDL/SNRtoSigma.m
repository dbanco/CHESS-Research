function sigma = SNRtoSigma(SNR,dataset)

% Compute noise variance needed for SNRs
[~,y_true,N,M,T,~,~] = sim_switch_multiscale_dl(0,dataset);
signal_power = mean(y_true(:).^2);
noise_power = signal_power/10^(SNR/10);
sigma = sqrt(noise_power);

end