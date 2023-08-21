% Test transposed operations
clear all
close all


%% Define parameters

% Length of intensity data (theta coordinate)
N = 55;
K = 10;
M = 50;
pad = M;
type = 'both';


% Define dictionary of Gaussian basis functions\
P.N = M;
P.K = K;   % Number of different basis functions 
P.sigmas = linspace(1/2,3,P.K); % Sigmas of basis functions

% Construct dictionary
A = zeros((N+M-1),(N+M-1)*K);
A0 = peakDictionary(P);

A0pad1 = padarray(A0,[N-1,0],0,'post');


A0pad2 = padarray(A0pad1,[2*pad 0],0,'post');
A0pad2pre = padarray(A0pad1,[2*pad 0],0,'pre');

A0ftpad1 = fft(A0pad1);
A0ftpad2 = fft(A0pad2);
A0ftpad2pre = fft(A0pad2pre);

figure
plot(A0pad1(:,5))
hold on
plot(A0(:,5))

for j = 1:K
    for i = 1:(N+M-1)
        ind1 = 1 + (i-1)*(N+M-1) + (j-1)*K*(N+M-1);
        ind2 =         i*(N+M-1) + (j-1)*K*(N+M-1);
        ind3 = i + (j-1)*(N+M-1);
        A(:,ind3) = circshift(A0pad1(:,j),i-1);
%         Apad = (:,ind3) = circshift(A0(:,j))
    end
end
x = zeros(N+M-1,K);
x(26,5) = 1;
x(66,5) = 1;

y1 = Ax_ft_1D(A0ftpad1,x);
y2 = A*x(:);

y3 = unpad(ifft(sum(bsxfun(@times,(A0ftpad1),fft(x)),2)),M-1,'pre');
y3pad = padarray(y3,M-1,0,'pre');
% y4pad = ifft(sum(bsxfun(@times,(A0ftpad2),fft(xpad)),2));
% y4 = unpad(y4pad,pad,type);
% y4pad2 = padarray(y3,[pad,0],0,type);

% norm(y1-y2)/norm(y2)
% norm(y1-y3)/norm(y1)
% norm(y2-y3)/norm(y2)
% norm(y4-y3)/norm(y3)
% norm(y4(61:100,:)-y3(61:100,:))/norm(y3(61:100,:))

% Aty1 = AtR_ft_1D(A0ftpad1,y1);
% Aty2 = reshape(A'*y1,[100,10]);
Aty3 = ifft(bsxfun(@times,conj(A0ftpad1),fft(y3pad)));
% Aty4pad = ifft(bsxfun(@times,conj(A0ftpad2),fft(y4pad2)));
% Aty4 = unpad(Aty4pad,pad,type);

% norm(Aty1(:)-Aty2(:))/norm(Aty2(:))
% norm(Aty1(:)-Aty3(:))/norm(Aty3(:))
% norm(Aty4(:)-Aty3(:))/norm(Aty3(:))
% norm(Aty4(61:100,:)-Aty3(61:100,:))/norm(Aty3(61:100,:))

figure(1)
subplot(2,1,1)
hold on
% plot(y1,'o')
% plot(y2)
plot(y3,'x-')
% plot(y4,'s')

subplot(2,1,2)
hold on
% plot(Aty1(:,10),'o')
% plot(Aty2(:,10))
plot(Aty3(:,5),'x-')
% plot(Aty4(:,10),'s')



