% SNR test script
% dataset = 'dissertation_adjust2';
dataset = 'pseudo-voigt_unmatched';
SNRs = [20,14,10,8,6];
SNR_estimates = zeros(numel(SNRs),1);
sigmas = zeros(numel(SNRs),1);

for i = 1:numel(sigmas)
    sigmas(i) = SNRtoSigma(SNRs(i),dataset);
    SNR_estimates(i) = SNRcheck(dataset,sigmas(i));
end

SNRs
SNR_estimates
sigmas