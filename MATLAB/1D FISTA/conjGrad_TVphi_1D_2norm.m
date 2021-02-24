function [Xk,cgIters] = conjGrad_TVphi_1D_2norm(A0ft_stack,B,X_init,YV,params)
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
N = size(A0ft_stack,1);

% Coefficeint Vectors
Xk = X_init;
Xk_lam = zeros(size(Xk));
for t = 1:size(Xk,3)
    Xk_lam(:,:,t) = Xk(:,:,t)*params.lambda1(t);
end

% Target Vectors
AtB = AtB_ft_1D_Time(A0ft_stack,B)./BnormSq2;

% Initial Residual
Rk = AtB - AtAx(A0ft_stack,Xk)./BnormSq2 +...
     + 2*params.lambda2*PtDtDPx(Xk) +...
     rho1*YV - rho1*Xk_lam;
Pk = Rk;

for i = 1:params.conjGradIter
    Apk = AtAx(A0ft_stack,Pk)./BnormSq2 + params.lambda2*PtDtDPx(Pk) + rho1*Pk;
    RkRk = sum(Rk(:).*Rk(:));
    alphak = RkRk/sum(Pk(:).*Apk(:));
    Xk = Xk + alphak*Pk;
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

function y = AtAx(A0ft_stack,X)
    y = AtB_ft_1D_Time(A0ft_stack,Ax_ft_1D_Time(A0ft_stack,X));
end

function y = PtDtDPx(X)
    N = size(X,1);
    y = PhiTranDiffTran_1D(DiffPhiX_1D(X),N);
end