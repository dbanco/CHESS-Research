function AtR = AtR_ft_3D( A0ft_stack, R )
%AtR_ft_3D Computes matrix vector product between A transposed and R
% Elementwise multiplies each basis function of A0ft_stack with fft2(R)
%
% Inputs:
% A0ft_stack - (n x m x l x t x r x w) array
% R - (n x m x l) array
% (n x m x l) is the size of the image and basis functions
% (t x r x w) indexes the basis function by theta variance and radial variance
%
% Outputs:
% AtR - (n x m x l x t x r x w) array


AtR = zeros(size(A0ft_stack));

R_ft = fftn(R);

for tv = 1:size(A0ft_stack,4)
    for rv = 1:size(A0ft_stack,5)
        for wv = 1:size(A0ft_stack,6)
            AtR(:,:,:,tv,rv,wv) = real(ifftn(A0ft_stack(:,:,:,tv,rv,wv).*R_ft));
        end
    end
end

end

