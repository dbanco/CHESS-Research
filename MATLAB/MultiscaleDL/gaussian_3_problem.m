function [y,D,X_true,N,M,T] = gaussian_3_problem
%% Construct 1D test problem Box and Gaussian
T = 60;
N = 128; M = N;
K = 1;
U = 3;
scaling = '2-norm';
d1 = gaussian_basis_wrap_1D(M,(M-1)/2,2,scaling);
d2 = gaussian_basis_wrap_1D(M,(M-1)/2,5,scaling);
d3 = gaussian_basis_wrap_1D(M,(M-1)/2,10,scaling);

D = zeros(1,M,3);
D(1,:,1) = d1;
D(1,:,2) = d2;
D(1,:,3) = d3;

% Ad1 = reSampleNu(N,d1,c1,c2,U);
% Ad2 = reSampleNu(N,d2,c1,c2,U);
% 
% AD = zeros(1,N,K*U);
% AD(1,:,1:U) = Ad1;
% AD(1,:,U+(1:U)) = Ad2;

a = 0.375*N;
Pos = round(a*sin( 2*pi*(1:T)/(T/2))+a);
% Pos2 = round(a*cos( 2*pi*(1:T)/(T/2))+a);

Wid = round((U-1)/2*sin( 2*pi*(1:T)/(T))+ (U-1)/2+1);
% Wid2 = round((U-1)/2*cos( 2*pi*(1:T)/(T))+ (U-1)/2+1); 

y = zeros(1,N,T);
X_true = zeros(1,N,K*U,T);
for t = 1:T
    y(:,:,t) = y(:,:,t) + circshift([D(:,:,Wid(t)), zeros(1,N-M)],Pos(t));
    X_true(1,Pos(t)+1,Wid(t),t) = 1;
end

end