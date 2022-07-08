function A0_stack = dictionaryFFT(P)
%dictionaryFFT Generates FFT of zero mean gaussian  
% basis function vectors with unit 2-norm
%
% Inputs:
% P.var_theta -vector of theta variances
% P.dtheta - difference in theta between each pixel
% P.num_theta - image size in theta direction 
%
% Outputs:
% A0ft_stack - FFT of dictionary atoms [N,K]

A0_stack = zeros(P.num_theta,numel(P.var_theta));
for t = 1:numel(P.var_theta)
    A0 = gaussian_basis_wrap_1D(P.num_theta, 1, P.var_theta(t),'2-norm');                    
    A0_stack(:,t) = fft(A0);
end

end

