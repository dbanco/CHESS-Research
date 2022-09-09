% Construct gaussian scale space
close all
K = 2;
P.Ny = 16; Ny = P.Ny;
P.Nx = 64; Nx = P.Nx;
P.stds = [0.8    0.8];
% P.stds = [0.8    0.8;
%           16    16];
P.means = [(P.Ny+1)/2 (P.Nx+1)/2]; 
D = dictionary2D(P);

%% Interpolator
a1 = 1.1;
a2 = 1.4;
j = 2;
S = 6;
Dpad = zeros(2*Ny,2*Nx,S);
psi = padarray(D(:,:,1),[Ny/2,Nx/2]);
x = linspace(-Nx,Nx,2*Nx);
y = linspace(-Ny,Ny,2*Ny);
z = psi;
[Y,X] = ndgrid(y,x);
F = griddedInterpolant(Y,X,z,'linear');
figure;
subplot(S,1,1)
imagesc(z)
title('Dilated D')
for j = 0:S
    z2 = F(a1^-j*Y,a2^-j*X);
    subplot(S+1,1,j+1)    
    Dpad(:,:,j+1) = circshift(z2,[-Ny/2,-Nx/2,0]);
    imagesc(Dpad(:,:,j+1))
end




%% Example sparse coding
X = zeros(Ny,Nx,S);
X(3,3,7) = 0;
X(10,10,4) = 0;
X(10,25,2) = 0;
X(10,1,3) = 0;
X(8,20,1) = 1;
X(5,3,5) = 0;
X(8,32,6) = 1;
Xpad = padarray(X,[Ny/2,Nx/2]); 

Xf = fft2(Xpad,2*Ny,2*Nx);
Df = fft2(Dpad,2*Ny,2*Nx);

Y = squeeze(sum(ifft2(Xf.*Df),3));

%% Generate full representation from just D1
% x = linspace(-Nx/2,Nx/2,Nx);
% y = linspace(-Ny/2,Ny/2,Ny);
x = linspace(-Nx,Nx,2*Nx);
y = linspace(-Ny,Ny,2*Ny);
% x = 1:2*Nx;
% y = 1:2*Ny;
[Yg,Xg] = ndgrid(y,x);

figure
subplot(S+1,2,1)
imagesc(Xpad(:,:,1))
title('X')
subplot(S+1,2,2)
imagesc(Xpad(:,:,1))
title('Dilated X')
Xhat = Xpad;

for j = 1:S
    i = j+1;
    z = Xpad(:,:,i);
    F = griddedInterpolant(Yg,Xg,z,'linear');
    z2 = F(a1^-j*Yg,a2^-j*Xg);
    subplot(S+1,2,2*i-1)
    imagesc(Xpad(:,:,i))
    subplot(S+1,2,2*i)
    imagesc(z2)
    Xhat(:,:,i) = z2;
end

Xhatf = fft2(Xhat,2*Ny,2*Nx);
Yhat = real(squeeze(sum(ifft2(Df(:,:,1).*Xhatf),3)));

exData = figure;
subplot(2,1,1)
imagesc(Y)
title('Dilated D Recon')
subplot(2,1,2)
imagesc(Yhat)
title('Dilated X Recon')

% Define phi 

