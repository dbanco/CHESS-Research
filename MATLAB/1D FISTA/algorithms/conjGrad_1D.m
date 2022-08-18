function [xk,cgIters] = conjGrad_1D(A0ft_stack,b,x_init,YV,params)
%conjGrad_TVx_1D Solves least squares
%
% Inputs:
%
% Outputs:
%


% ADMM penalty parameter
rho1 = params.rho1;

% Coefficeint Vectors
xk = x_init;

% Target Vectors
AtB = AtR_ft_1D(A0ft_stack,b);

% Initial Residual
Rk = AtB - AtR_ft_1D(A0ft_stack,Ax_ft_1D(A0ft_stack,xk)) +...
     rho1*YV  - rho1*xk;
Pk = Rk;

for i = 1:params.conjGradIter
    Apk = AtR_ft_1D(A0ft_stack,Ax_ft_1D(A0ft_stack,Pk)) + rho1*Pk;
    RkRk = sum(Rk(:).*Rk(:));
    alphak = RkRk/sum(Pk(:).*Apk(:));
    xk = xk + alphak*Pk;
    Rkp1 = Rk - alphak*Apk;
    if norm(Rkp1(:)) < params.cgEpsilon
        break;
    end
    betak = sum(Rkp1(:).*Rkp1(:))/RkRk;
    Pk = Rkp1 + betak*Pk;
    Rk = Rkp1;
end
cgIters = i;
end
