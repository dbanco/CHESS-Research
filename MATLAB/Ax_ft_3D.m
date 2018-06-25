function Ax = Ax_ft_3D( A0ft_stack, x )
%Ax_ft_2D Computes matrix vector product between A and x
% Elementwise multiplies each basis function of A0ft_stack with fft2(x)
%
% Inputs:
% A0ft_stack - (n x m x l x t x r x w) array
% x - (n x m x l x t x r x w ) array
% (n x m x l) is the size of the image and basis functions
% (t x r x w) indexes the basis function by theta variance and radial variance
%
% Outputs:
% Ax - (n x m x l) array

Ax = zeros(size(A0ft_stack,1),size(A0ft_stack,2),size(A0ft_stack,3));

x_ft = fftn(x);

for tv = 1:size(A0ft_stack,4)
    for rv = 1:size(A0ft_stack,5)
        for wv = 1:size(A0ft_stack,6)
            Ax = Ax + real(ifftn(A0ft_stack(:,:,:,tv,rv,wv).*x_ft(:,:,:,tv,rv,wv)));
        end
    end
end

end

