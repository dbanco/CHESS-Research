function [y,ydu,X_true,N,T] = upDwn_translate_problem
%% Construct 1D test problem Box and Gaussian
T = 30;
N = 128;
f = zeros(1,N);
x = 1:N-16;
f(9:N-8) = sinc(8*x/N - pi) + 0.25;
plot(f)
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

yd = decimate(reshape(Pnrm(f),[1,N,1,1]),5);
ff = Pnrm(yd(:,1:N/4,3));

ydu = upDwnSample(ff,2);

plotDictionary(ydu)

Pos = round(48*sin( 2*pi*(1:T)/12)+48);
% plot(Pos)

Wid = round(2*sin( 2*pi*(1:T)/15)+3);
% plot(Wid)

y = zeros(1,N,T);
X_true = zeros(1,N,5,T);
for t = 1:T
    y(:,:,t) = y(:,:,t) + circshift(ydu(:,:,Wid(t)),Pos(t));
    X_true(1,Pos(t)+1,Wid(t),t) = 1;
end



end