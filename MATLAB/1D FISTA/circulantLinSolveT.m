function [x_hat] = circulantLinSolveT( A0ft_stack,B,Yk,Vk,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

Bnormsq = sum(B(:).^2);
B_ft = fft(B);

rho = params.rho*Bnormsq;

r_ft = (A0ft_stack'.').*repmat(b_ft,[1,20]) + rho.*fft(yk-vk);

coef = diag(A0ft_stack*r_ft.')./...
      (rho + diag(A0ft_stack*A0ft_stack'));
  
x_ft = (r_ft - repmat(coef,[1,20]).*(A0ft_stack'.') )/(rho);

x_hat = real(ifft(x_ft));
end

