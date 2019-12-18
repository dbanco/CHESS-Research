function Ax = Ax_ft_2D( A0ft_stack, x )
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

Ax = zeros(size(A0ft_stack,1),size(A0ft_stack,2));

x_ft = fft2(x);

for tv = 1:size(A0ft_stack,3)
    for rv = 1:size(A0ft_stack,4)
        Ax = Ax + real(ifft2(A0ft_stack(:,:,tv,rv).*x_ft(:,:,tv,rv)));
    end
end


end

