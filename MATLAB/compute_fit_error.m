function rel_fit_error = compute_fit_error( A0ft_stack,x,b )
%compute_fit_error Summary of this function goes here
%   Detailed explanation goes here
fit = Ax_ft_2D(A0ft_stack,x);
rel_fit_error = norm( b(:)-fit(:) )/norm(b(:));
end

