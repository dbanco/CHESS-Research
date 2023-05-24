function Axf = Ax_gpu( Af, xf )
%Ax_ft_2D Computes matrix vector product between A and x
% Elementwise multiplies each basis function of Af with xf
%
% Inputs:
% Af - (N x M x K x T) array
% Xf - (N x M x K x T) array
% (N x M) is the size of the image and basis functions
%
% Outputs:
% Axf - (N x M x 1 x T) array

Axf = sum(pagefun(@times, Af, xf),3);


end

