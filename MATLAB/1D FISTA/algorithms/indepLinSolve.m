function [y1kp1] = indepLinSolve(A0ft_stack,b,xkp1,v1k,params )
%circulantLinSolve Solves M independent linear systems using
%Sherman-Morrison formula as described in "Efficient Convolutional Sparse
%Coding" Brendt Wohlberg

rho = params.rho1;
zMask = params.zeroMask;
Wb = forceMaskToZero(b,zMask);

r = Wb + rho*(Ax_ft_1D(A0ft_stack,xkp1) + v1k);
coef = forceMaskToZero(ones(size(v1k)),zMask) + rho*ones(size(v1k)) ;
y1kp1 = linsolve(diag(coef),r);

if params.plotProgress
    figure(3)
    hold off
    fit = Ax_ft_1D(A0ft_stack,y1kp1);
    plot(b)
    hold on
    plot(fit)
    legend('b','fit')
end

end

