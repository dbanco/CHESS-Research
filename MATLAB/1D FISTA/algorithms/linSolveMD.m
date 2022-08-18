function [x_hat] = linSolveMD(b,y1k,v1k,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

bnormsq = sum(b(:).^2);
rho = params.rho1*bnormsq;
Wb = forceMaskToZero(b,params.zMask);

r = Wb + rho*(Ax_ft_1D(A0ft_stack,x) + v1k;
coef = ( diag(forceMaskToZero(ones(size(y1k)))) + rho*eye(numel(y1k)) )*y1k;
x_hat = r./coef;


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

