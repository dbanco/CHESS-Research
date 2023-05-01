% Test transposed operations
clear all
close all


%% Define parameters

% Length of intensity data (theta coordinate)
N = 100;
K = 10;

% Define dictionary of Gaussian basis functions\
P.N = N;
P.K = K;   % Number of different basis functions 
P.sigmas = linspace(1/2,16,P.K); % Sigmas of basis functions

% Construct dictionary
A = zeros(P.N,P.N*P.K);
% Construct dictionary
A0ft = peakDictionaryFFT(P);
A0 = peakDictionary(P);

for j = 1:K
    for i = 1:P.N
        ind1 = 1 + (i-1)*N + (j-1)*K*N;
        ind2 =         i*N + (j-1)*K*N;
        ind3 = i + (j-1)*N;
        A(:,ind3) = circshift(A0(:,j),i-1);
    end
end
x = randn(P.N,P.K);

y1 = Ax_ft_1D(A0ft,x);
y2 = A*x(:);
y3 = ifft(sum(bsxfun(@times,(A0ft),fft(x)),2));

norm(y1-y2)/norm(y2)
norm(y1-y3)/norm(y1)
norm(y2-y3)/norm(y2)

Aty1 = AtR_ft_1D(A0ft,y1);
Aty2 = A'*y1;
Aty3 = ifft(bsxfun(@times,conj(A0ft),fft(y3)));

norm(Aty1(:)-Aty2)/norm(Aty2)
norm(Aty1(:)-Aty3(:))/norm(Aty3(:))

figure(1)
subplot(2,1,1)
hold on
plot(y1,'o')
plot(y2)
plot(y3,'x')

subplot(2,1,2)
hold on
plot(Aty1(:),'o')
plot(Aty2(:))
plot(Aty3(:),'x')