function Ax = Ax_ft_2D( Af, xf )
%Ax_ft_2D Computes matrix vector product between A and x
% Elementwise multiplies each basis function of A0ft_stack with fft2(x)
%
% Inputs:
% A0ft_stack - (n x m x t x r) array
% x - (n x m x t x r) array
% (n x m) is the size of the image and basis functions
% (t x r) indexes the basis function by theta variance and radial variance
%
% Outputs:
% Ax - (n x m) array

Ax = zeros(size(Af,1),size(A0ft_stack,2));

if 
Axf = bsxfun(@times, conj(Af), xf);


end

