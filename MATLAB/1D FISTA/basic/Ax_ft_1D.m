function Ax = Ax_ft_1D( A0ft_stack, x )
%Ax_ft_2D Computes matrix vector product between A and x
% Elementwise multiplies each basis function of A0ft_stack with fft2(x)
%
% Inputs:
% A0ft_stack - (n x t) array
% x - (n x t) array
% (n x 1) is the size of the image and basis functions
% (t x 1) indexes the basis function by theta variance
%
% Outputs:
% Ax - (n x m) array

Ax = zeros(size(A0ft_stack,1),1);

x_ft = fft(x);

for tv = 1:size(x,2)
        Ax = Ax + real(ifft(A0ft_stack(:,tv).*x_ft(:,tv)));
end

end

