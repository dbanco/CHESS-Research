function [yn,y,N,M,T] = generate_signal_pair(sigma)
%% Construct 1D test problem Gaussian and linear 
T = 30;
N = 55; 
M = 55;

if nargin < 1
    sigma = 0;
%     sigma = 0.015;
end

Pos = linspace(8,40,T);
Pos2 = linspace(45,10,T);

sig = linspace(3,8,T);
width = round(linspace(20,5,T));

% Position in time
y = zeros(N,T);

center = (N-1)/2 + 1;

for t = 1:T
    tooth = zeros(N,1);
    if mod(width(t),2)
        tooth(center-floor(width(t)/2): center+floor(width(t)/2)) = linspace(0,1,width(t))/norm(linspace(0,1,width(t)));
    else
        tooth(center-floor(width(t)/2): center+floor(width(t)/2)-1) = linspace(0,1,width(t))/norm(linspace(0,1,width(t)));
    end
    y(:,t) = gaussian_basis_wrap_1D(N,Pos(t),sig(t),'2-norm') + circshift(tooth,round(Pos2(t)-center)) ;
end

yn = y + randn(N,T)*sigma;
% yn(yn<0)=0;

% snr = norm(y(:))/norm(y(:)-yn(:))
% Reduce data to a time subset
% trange = 1:60;
% yn = yn(:,:,trange);
% y = y(:,:,trange);
% Xtrue = Xtrue(:,:,:,trange);
% T = numel(trange);

% figure
% imagesc(squeeze(yn))
% 
% figure
% plot(y(:,1))

end