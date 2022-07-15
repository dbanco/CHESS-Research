function [y,ydu,X_true,N,T,scales] = upDwn_double_multirate_problem
%% Construct 1D test problem Box and Gaussian
T = 60;
N = 128;
K = 2;
scales = [1 2 3 3;
          3 3 2 1];
U = size(scales,2);
ydu1 = multiRateSinc(N,scales);
ydu2 = multiRateShape(N,scales);

% plotDictionary(ydu2)
Pos = round(48*sin( 2*pi*(1:T)/(T/2))+48);
Pos2 = round(48*cos( 2*pi*(1:T)/(T/2))+48);
% plot(Pos)

% Should oscillate between 1 and U+1
Wid = round((U-1)/2*sin( 2*pi*(1:T)/(T))+ (U-1)/2+1);
Wid2 = round((U-1)/2*cos( 2*pi*(1:T)/(T))+ (U-1)/2+1); 
% plot(Wid)

y = zeros(1,N,T);
X_true = zeros(1,N,K*U,T);
for t = 1:T
    y(:,:,t) = y(:,:,t) + circshift(ydu1(:,:,Wid(t)),Pos(t));
    X_true(1,Pos(t)+1,Wid(t),t) = 1;
    y(:,:,t) = y(:,:,t) + circshift(ydu2(:,:,Wid2(t)),Pos2(t));
    X_true(1,Pos2(t)+1,Wid2(t)+U,t) = 1;
end

ydu = zeros(1,N,2*U);
ydu(:,:,1:U) = ydu1;
ydu(:,:,U+1:end) = ydu2;

end