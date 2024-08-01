function [yn,y,N,M,T] = gaus_example_matched_multiscale_dl()
%% Construct 1D test problem Gaussian and linear 
T = 50;
N = 105; M = 55;

if nargin < 1
    sigma = 0.01;
end

% Position in time
position = 15+[linspace(20,50,20),30+20.*linspace(1,0.1,30).*cos((0:29)/3)];
width = (70-fliplr(position))/3;

max_width = max(width);

% Model Setup
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
KJ = K*J;
center = (M+1)/2;

Dtrue = zeros(1,M,K);
Dtrue(1,:,1) = gaussian_basis_wrap_1D(M,center,max_width,'2-norm');
ADtrue = padarray(reSampleCustomArrayCenter(N,Dtrue,scales,center),[0,M-1,0,0],0,'post');
ADf = fft2(ADtrue);

position = round(position)+center;
width = round(width);

Xtrue = zeros(1,N+M-1,KJ,T);
for t = 1:T
    Xtrue(1,position(t),width(t),t) = 1;
end

Xf = fft2(Xtrue);
y = squeeze(unpad(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'),M-1,'pre'));

yn = y + randn(N,T)*sigma;

figure(2)
imagesc(yn)

end