function [y,yu,X_true,N,T,U] = upSample_translate_problem
%% Construct 1D test problem Box and Gaussian
T = 30;
N = 128;
f = zeros(1,N);
x = 1:N-16;
f(9:N-8) = sinc(8*x/N - pi) + 0.25;
plot(f)
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

U = 3;
yd = decimate(reshape(Pnrm(f),[1,N,1,1]),U);
Nd = N/(2^(U-1));
ff = reshape(yd(:,1:Nd,U),[1,Nd,1,1]);
yu = upSample(Pnrm(ff),U);
plotDictionary(yu)
norm(yd(:,:,3))
norm(yd(:,:,2))

Pos = round(N*3/8.*sin( 2*pi*(1:T)/12)+N*3/8);
% plot(Pos)

Wid = round(floor(U/2)*sin( 2*pi*(1:T)/15)+ceil(U/2));
% plot(Wid)

y = zeros(1,N,T);
X_true = zeros(1,N,U,T);
for t = 1:T
    y(:,:,t) = y(:,:,t) + circshift(yu(:,:,Wid(t)),Pos(t));
    X_true(1,Pos(t)+1,Wid(t),t) = 1;
end



end