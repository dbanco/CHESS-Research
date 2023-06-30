function [objOF, objHS, sys] = HSobjectivePaper(Fx,Fy,Ft,u,v,K,smoothness)

KU = size(Fx,2);
U = KU/K;

[nablaU2,nablaV2] = laplaceOp(u,v,K,U);

objOF = sum((Fy.*v + Fx.*u + Ft).^2,'all');
objHS = sum((nablaV2).^2,'all') + sum((nablaU2).^2,'all');

b = [-Fx(:).*Ft(:); -Fy(:).*Ft(:)];

FxFy = Fx(:).*Fy(:);
Au = Fx(:).^2.*u(:) + FxFy.*v(:) - smoothness*nablaU2(:) ;
Av = Fy(:).^2.*v(:) + FxFy.*u(:) - smoothness*nablaV2(:) ;

Ax = [Au; Av];

sys = b-Ax;




end
