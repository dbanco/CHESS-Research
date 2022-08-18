function [Xk,cgIters] = conjGrad_TVphi_1D(A0ft_stack,B,Bnorms,X_init,YV,ZU,params,zMask)
%conjGrad_TVx_1D Solves least squares
%
% Inputs:
%
% Outputs:
%

% ADMM penalty parameter
rho1 = params.rho1;
rho2 = params.rho2;
N = size(A0ft_stack,1);

% Coefficeint Vectors
Xk = X_init;

% Target Vectors
AtB = AtB_ft_1D_Time(A0ft_stack,B,Bnorms);
PtDtZ = PhiTranDiffTran_1D(ZU,N);

% Initial Residual
Rk = AtB - AtAx(A0ft_stack,Xk,Bnorms) +...
     rho2*PtDtZ - rho2*PtDtDPx(Xk) +...
     rho1*YV - rho1*Xk;
Pk = Rk;

for i = 1:params.conjGradIter
    Apk = AtAx(A0ft_stack,Pk,Bnorms) + rho2*PtDtDPx(Pk) + rho1*Pk;
    RkRk = sum(Rk(:).*Rk(:));
    alphak = RkRk/sum(Pk(:).*Apk(:));
    Xk = Xk + alphak*Pk;
    Xk = forceMaskToZeroArray(Xk,zMask);
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

function y = AtAx(A0ft_stack,X,Bnorms)
    y = AtB_ft_1D_Time(A0ft_stack,Ax_ft_1D_Time(A0ft_stack,X,Bnorms),Bnorms);
end

function y = PtDtDPx(X)
    N = size(X,1);
    y = PhiTranDiffTran_1D(DiffPhiX_1D(X),N);
end