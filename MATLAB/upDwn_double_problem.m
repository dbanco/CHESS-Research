function [y,ydu,X_true,N,T] = upDwn_double_problem
%% Construct 1D test problem Box and Gaussian
T = 60;
N = 128;
U = 2;
K = 2;
ydu1 = multiScaleSinc(N,U+1);
ydu2 = multiScaleShape(N,U+1);

% plotDictionary(ydu2)

Pos = round(48*sin( 2*pi*(1:T)/(T/2))+48);
Pos2 = round(48*cos( 2*pi*(1:T)/(T/2))+48);
% plot(Pos)

Wid = round(U*sin( 2*pi*(1:T)/(T))+ U+1);
Wid2 = round(U*cos( 2*pi*(1:T)/(T))+ U+1); 
% plot(Wid)

y = zeros(1,N,T);
X_true = zeros(1,N,K*(2*U+1),T);
for t = 1:T
    y(:,:,t) = y(:,:,t) + circshift(ydu1(:,:,Wid(t)),Pos(t));
    X_true(1,Pos(t)+1,Wid(t),t) = 1;
    y(:,:,t) = y(:,:,t) + circshift(ydu2(:,:,Wid2(t)),Pos2(t));
    X_true(1,Pos2(t)+1,Wid2(t)+(2*U+1),t) = 1;
end

ydu = zeros(1,N,2*(2*U+1));
ydu(:,:,1:(2*U+1)) = ydu1;
ydu(:,:,(2*U+1)+1:end) = ydu2;

end