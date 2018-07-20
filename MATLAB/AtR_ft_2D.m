function AtR = AtR_ft_2D( A0ft_stack, R )
%AtR_ft_2D Computes matrix vector product between A transposed and R
% Elementwise multiplies each basis function of A0ft_stack with fft2(R)
%
% Inputs:
% A0ft_stack - (n x m x t x r) array
% R - (n x m) array
% (n x m) is the size of the image and basis functions
% (t x r) indexes the basis function by theta variance and radial variance
%
% Outputs:
% AtR - (n x m x t x r) array


AtR = zeros(size(A0ft_stack));

R_ft = fft2(R);

for tv = 1:size(A0ft_stack,3)
    for rv = 1:size(A0ft_stack,4)
         y = real(ifft2(A0ft_stack(:,:,tv,rv).*R_ft));
        AtR(:,:,tv,rv) = shift2D(y,2,2);
    end
end

end

