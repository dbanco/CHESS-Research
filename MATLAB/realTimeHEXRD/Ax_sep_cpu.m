function Axf = Ax_sep_cpu( Af1, Af2, xf )
%Ax_ft_2D Computes matrix vector product between A and x
% Elementwise multiplies each basis function of Af with xf
%
% Inputs:
% Af1 - (N x 1 x K1 x 1  x T) array
% Af2 - (1 x M x 1  x K2 x T) array
% Xf -  (N x M x K1 x K2 x T) array
% (N x M) is the size of the image and basis functions
%
% Outputs:
% Axf - (N x M x 1 x T) array

Axf = zeros(size(xf));

for i = 1:size(xf,3)
    conv1 = bsxfun(@times, Af1(:,:,i,:), xf(:,:,i,:));
    for j = 1:size(xf,4)
        Axf(:,:,i,j) = bsxfun(@times,Af2(:,:,:,j),conv1(:,:,:,j));
    end
end
Axf = sum(Axf,[3,4])

end

