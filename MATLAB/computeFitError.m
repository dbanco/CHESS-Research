function relFitError = computeFitError( A0ft_stack,x,b )
%computeFitError
fit = Ax_ft_2D(A0ft_stack,x);
relFitError = norm( b(:)-fit(:) )/norm(b(:));
end

