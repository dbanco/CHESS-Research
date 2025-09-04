function SNR = SNRcheck(dataset,sigma_val)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[y,y_true,N,M,T] = sim_switch_multiscale_dl(sigma_val,dataset);
signal_power = mean(y_true(:).^2);
noise_power = mean((y(:)-y_true(:)).^2);
SNR = 10*log10(signal_power/noise_power);
noiseTheor = sigma_val^2;


end