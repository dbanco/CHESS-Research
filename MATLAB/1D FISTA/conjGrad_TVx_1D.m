function Xk = conjGrad_TVx_1D(A0ft_stack,B,X_init,YV,ZU,params)
%conjGrad_TVx_1D Solves least squares
%
% Inputs:
%
% Outputs:
%

% Data normalizing
BnormSq1 = (sum(B.^2,1));
BnormSq2 = reshape(BnormSq1,[1,1,numel(BnormSq1)]);

% ADMM penalty parameter
rho1 = params.rho1;
rho2 = params.rho2;

% Coefficeint Vectors
Xk = X_init;

% Target Vectors
AtB = AtB_ft_1D_Time(A0ft_stack,B)./BnormSq2;
DtZ = DiffTranB_1D(ZU);
DtDx(Xk);

% Initial Residual
Rk = AtB - AtAx(A0ft_stack,Xk)./BnormSq2 +...
     rho2*DtZ - rho2*DtDx(Xk) +...
     rho1*YV  - rho1*Xk;
Pk = Rk;

for i = 1:params.conjGradIter
    Apk = AtAx(A0ft_stack,Pk)./BnormSq2 + rho2*DtDx(Pk) + rho1*Pk;
    RkRk = sum(Rk(:).*Rk(:));
    alphak = RkRk/sum(Pk(:).*Apk(:));
    Xk = Xk + alphak*Pk;
    Rkp1 = Rk - alphak*Apk;
    betak = sum(Rkp1(:).*Rkp1(:))/RkRk;
    Pk = Rkp1 + betak*Pk;
    Rk = Rkp1;
end

end

function y = AtAx(A0ft_stack,X)
    y = AtB_ft_1D_Time(A0ft_stack,Ax_ft_1D_Time(A0ft_stack,X));
end

function y = DtDx(X)
    y = DiffTranB_1D(DiffX_1D(X));
end