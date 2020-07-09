function AtR = AtR_ft_1D( A0ft_stack, R )
%AtR_ft_1D Computes matrix vector product between A transposed and R
% Elementwise multiplies each basis function of A0ft_stack with fft(R)
%
% Inputs:
% A0ft_stack - (n x t) array
% R - (n x 1) array
% (n x 1) is the size of the image and basis functions
% (t x 1) indexes the basis function by theta variance 
%
% Outputs:
% AtR - (n x m x t x r) array


AtR = zeros(size(A0ft_stack));

R_ft = fft(R);

for tv = 1:size(A0ft_stack,2)
        y = real(ifft(A0ft_stack(:,tv).*R_ft(:)));
        AtR(:,tv) = y;
end

AtR = shift1D_array(AtR,2);

end

