function [y,AD,D,X_true,N,M,T,scales,c1,c2] = upDwn_double_multirate_problem
%% Construct 1D test problem Box and Gaussian
T = 60;
N = 128; M = 32;
K = 2;
c1 = 5;
c2 = 3;
scales = [ c2 1 c1;
           c1 1 c2];
U = size(scales,2);
d1 = multiRateSinc(M);
d2 = multiRateShape(M);

D = zeros(1,M,K);
D(1,:,1) = d1;
D(1,:,2) = d2;

Ad1 = reSampleNu(N,d1,c1,c2,U);
Ad2 = reSampleNu(N,d2,c1,c2,U);

AD = zeros(1,N,K*U);
AD(1,:,1:U) = Ad1;
AD(1,:,U+(1:U)) = Ad2;

% figure(1)
% plot(d1)
% hold on
% plot(d2)

%%

a = 0.375*N;
Pos = round(a*sin( 2*pi*(1:T)/(T/2))+a);
Pos2 = round(a*cos( 2*pi*(1:T)/(T/2))+a);

Wid = round((U-1)/2*sin( 2*pi*(1:T)/(T))+ (U-1)/2+1);
Wid2 = round((U-1)/2*cos( 2*pi*(1:T)/(T))+ (U-1)/2+1); 

y = zeros(1,N,T);
X_true = zeros(1,N,K*U,T);
for t = 1:T
    y(:,:,t) = y(:,:,t) + circshift(Ad1(:,:,Wid(t)),Pos(t));
    X_true(1,Pos(t)+1,Wid(t),t) = 1;
    y(:,:,t) = y(:,:,t) + circshift(Ad2(:,:,Wid2(t)),Pos2(t));
    X_true(1,Pos2(t)+1,Wid2(t)+U,t) = 1;
end

end