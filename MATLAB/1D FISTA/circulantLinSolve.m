function [x_hat] = circulantLinSolve( A0ft_stack,b,yk,vk,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

% bnormsq = sum(b(:).^2);
b_ft = fft(b);
K = size(A0ft_stack,2);
rho = params.rho1;

r_ft = (A0ft_stack'.').*repmat(b_ft,[1,K]) + rho*fft(yk-vk);
coef = diag(A0ft_stack*r_ft.')./...
      (rho + diag(A0ft_stack*A0ft_stack'));
x_ft = (r_ft - repmat(coef,[1,K]).*(A0ft_stack'.') )/rho;
x_hat = real(ifft(x_ft));

if params.plotProgress
    figure(2)
    hold off
    fit = Ax_ft_1D(A0ft_stack,x_hat);
    plot(b)
    hold on
    plot(fit)
    legend('b','fit')
end

end

