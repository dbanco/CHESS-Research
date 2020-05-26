function [x_hat] = circulantLinSolve( A0ft_stack,b,yk,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg
%

rho = params.rho;

b_ft = fft(b);
y_ft = fft(yk);

r_ft = (A0ft_stack'.').*repmat(b_ft,[1,20]) + rho*y_ft;
coef = diag(A0ft_stack*r_ft.')./(rho + diag(A0ft_stack*A0ft_stack'));
x_ft = (r_ft - repmat(coef,[1,20]).*(A0ft_stack'.'))/rho;
x_hat = real(ifft(x_ft));
end

