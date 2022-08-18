function [x_hat] = circulantLinSolveMD( A0ft_stack,b,y0k,v0k,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

b_ft = fft(b);
K = size(A0ft_stack,2);


r_ft = (A0ft_stack'.').*repmat(b_ft,[1,K]) + fft(y0k-v0k);
coef = diag(A0ft_stack*r_ft.')./...
      (1 + diag(A0ft_stack*A0ft_stack'));
x_ft = (r_ft - repmat(coef,[1,K]).*(A0ft_stack'.') );
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

