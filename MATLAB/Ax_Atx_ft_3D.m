function y = Ax_Atx_ft_3D( x, type, A0ft_stack )
%Ax_Atx_ft_3D Summary of this function goes here
%   Detailed explanation goes here

if strcmp(type,'notransp')
    Ax = zeros(size(A0ft_stack,1),size(A0ft_stack,2));
    x_ft = fft2(reshape(x,size(A0ft_stack)));
    
    for tv = 1:size(A0ft_stack,3)
        for rv = 1:size(A0ft_stack,4)
            Ax = Ax + real(ifft2(A0ft_stack(:,:,tv,rv).*x_ft(:,:,tv,rv)));
        end
    end
    y = reshape(Ax,numel(Ax),1);
    
elseif strcmp(type,'transp')
    AtR = zeros(size(A0ft_stack));
    R_ft = fft2(reshape(x,size(AtR,1),size(AtR,2)));

    for tv = 1:size(A0ft_stack,3)
        for rv = 1:size(A0ft_stack,4)
            AtR(:,:,tv,rv) = real(ifft2(A0ft_stack(:,:,tv,rv).*R_ft));
        end
    end
    y = reshape(AtR,numel(AtR),1);
end

end

