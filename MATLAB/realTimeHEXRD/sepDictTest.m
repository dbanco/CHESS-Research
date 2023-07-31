% Separable dictionary test script
P.N1 = 30;
P.N2 = 300;
P.mu1 = 15;
P.mu2 = 150;
P.K1 = 10;
P.K2 = 15;
P.K = P.K1*P.K2;
P.sigma1 = linspace(0.5,  4,    P.K1);
P.sigma2 = linspace(0.5,  12,   P.K2); % was 12
P.basis = 'norm2';

D = dictionary2D(P);
[sepD1,sepD2] = sepDictionary2D(P);

X = zeros(P.N1,P.N2,P.K1*P.K2);
% X(4,100,11) = 1;
% X(9,100,2) = 1;
X(1,1,4) = 1;
% X(15,10,122) = 1;
% X(20,100,140) = 1;

X2 = reshape(X,[P.N1,P.N2,P.K1,P.K2]);
X1f = fft2(X);
X2f = fft2(X2);

tic
y1 = real(ifft2(Ax_cpu(fft2(D),X1f)));
toc

tic
y2 = real(ifft2(Ax_sep_cpu(fft2(sepD1),fft2(sepD2),X2f)));
toc

figure(3)
imagesc(y1)
figure(4)
imagesc(y2)

figure(1)
imagesc(D(:,:,1))
figure(2)
imagesc(sepD1(:,:,1,:)*sepD2(:,:,:,1))

norm( D(:,:,1) - sepD1(:,:,1,:)*sepD2(:,:,:,1) )

norm(y1-y2)