function Atxf = Atx_cpu( Af, xf )
%Ax_ft_2D Computes matrix vector product between A and x
% Elementwise multiplies each basis function of Af with xf
%
% Inputs:
% Af - (N x M x K x T) array
% Xf - (N x M x K x T) array
% (N x M) is the size of the image and basis functions
%
% Outputs:
% Atxf - (N x M x 1 x T) array

Atxf = sum(bsxfun(@times, conj(Af), xf),4);


end

