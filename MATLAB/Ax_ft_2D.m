function Ax = Ax_ft_2D( A0ft_stack, x )
%AX_ft_2D Summary of this function goes here
%   Detailed explanation goes here

Ax = zeros(size(A0ft_stack,1),size(A0ft_stack,2));

x_ft = fft2(x);

for tv = 1:size(A0ft_stack,3)
    for rv = 1:size(A0ft_stack,4)
        Ax = Ax + real(ifft2(A0ft_stack(:,:,tv,rv).*x_ft(:,:,tv,rv)));
    end
end

end

