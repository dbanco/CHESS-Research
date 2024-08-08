sigmas = 0.01:0.01:0.05;
lambda_sel = [18,20,22,20,19];

% Load outputs
for j = 1:5 
    load(['C:\Users\dpqb1\Documents\Outputs\toy1_exp1_sig_',...
           num2str(j),'\output_',num2str(lambda_sel(j)),'.mat'])
    
    % Compute SNR
    signal_power = squeeze(sum(outputs.y.^2,2))./outputs.N;
    noise_power = sigmas(j)^2;
    
    snr = mean(signal_power./noise_power)
end