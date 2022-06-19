% LSQR test and matrix operations test
U = 3; N = 64; K = 2; T = 20;
y = zeros(1,N,1,T); y(1,N/4:N/2,1,1) = 1;
x = zeros(1,N,K*U,T); x(1,16,2,1) = 1; x(1,32,3,1) = 0;
d = 0.5*ones(1,N,K);

% % Matrix operations
% % Ax operation (A is decimation matrix)
% Ad = decimate(d,U);
% % XAd operation
% XAd = ifft2(sum(bsxfun(@times,fft2(x),fft2(Ad)),3),'symmetric');
% figure(1)
% plot(XAd(1,:,1,1))
% 
% % Transposed matrix operations
% % Xty operation
% Xty = ifft2(sum(bsxfun(@times, fft2(x), fft2(y)), 4),'symmetric');
% % AtXty operation
% AtXty = decimateTrans(Xty,U);
% figure(2)
% imagesc(squeeze(Xty))
% figure(3)
% imagesc(squeeze(AtXty))

%% LSQR
A = @(in,type) AfuncDecim(in,type,U,x,N,T,K); 
tol = 1e-6;
maxIters = 200;
lsqrSOL(N*T,N*K,A,y(:),0,tol,tol,tol,maxIters,1)