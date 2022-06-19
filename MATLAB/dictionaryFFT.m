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

K = numel(P.variances);
A0_stack = zeros(P.N,K);
for i = 1:K
    A0 = gaussian_basis_wrap_1D(P.N, 1, P.variances(i),'2-norm');                    
    A0_stack(:,i) = fft(A0);
end

end

