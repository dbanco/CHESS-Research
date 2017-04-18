function AtR = AtR_ft_2D( A0ft_stack, R )
%AX_ft_2D Summary of this function goes here
%   Detailed explanation goes here

AtR = zeros(size(A0ft_stack));

R_ft = fft2(R);

for tv = 1:size(A0ft_stack,3)
    for rv = 1:size(A0ft_stack,4)
        AtR(:,:,tv,rv) = real(ifft2(A0ft_stack(:,:,tv,rv).*R_ft));
    end
end

end

