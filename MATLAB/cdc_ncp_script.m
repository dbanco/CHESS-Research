% Dictionary Constrained Nonnegative CP

% Load MMPAD data
data_dir = 'D:\MMPAD_data_nr1\ring1_zero';
t = 1:4:200;
T = numel(t);
for i = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t(i)),'.mat']))
    if t == 1
        [N,M] = size(polar_image);
        S0 = zeros(N,M,T);
    end
    S0(:,:,i) = polar_image/norm(polar_image,'fro');
end

% Set up cbpdndl parameters
lambda = 0.01;
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 50;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(S0,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;

param.R = 20;
param.maxIters = 5;
param.ncp_iters = 50;       
param.lambda = lambda;
param.gamma = 0.01;
param.stepSize = 0.1;
param.verbose = 1;

% Initial 2theta dictionary
K1 = 9;
opt1 = opt;
opt1.DictFilterSizes = [ones(1,K1);
                       8,8,8,16,16,16,32,32,32];
P1.N = max(opt1.DictFilterSizes(:));
P1.variances = logspace(log10(0.5),log10(P1.N),K1).^2;
P1.mean = 4;
P1.basis = 'norm2';
D1 = zeros(1,P1.N,K1);
for i = 1:K1
    D1(1,:,:) = dictionary(P1);
end

% Initial eta dictionary
K2 = 13;
opt2 = opt;
opt2.DictFilterSizes = [ones(1,K2);
                       16,16,16,32,32,32,64,64,64,130,130,261,261];
P2.N = max(opt2.DictFilterSizes(:));
P2.variances = logspace(log10(0.5),log10(P2.N),K2).^2;
P2.mean = 8;
P2.basis = 'norm2';
D2 = zeros(1,P2.N,K2);
for i = 1:K2
    D2(1,:,:) = dictionary(P2);
end

%% Solve cdc_ncp
[CP, D1, D2, X1, X2] = cdc_ncp(S0,D1,D2,param,opt1,opt2);

%% View result
figure(1)
n = 1;
F = double(CP);
subplot(2,1,1)
imagesc(F(:,:,n))
title('recon')

subplot(2,1,2)
imagesc(S0(:,:,n))
title('data')

%% Visualize factors in time order
% Sort by time occurence
U1 = CP.U{1};
U2 = CP.U{2};
U3 = CP.U{3};

% Compute time center of mass
[T,R]= size(U3);
timeW = 1:T;
timeCoM = zeros(R,1);
for i = 1:R
    timeCoM(i) = timeW*U3(:,i)/sum(U3(:,i),'all');
end
[~,ind] = sort(timeCoM);

figure(2) % View angle components
for i = 1:R
    subplot(5,4,i)
    j = ind(i);
    angComp = CP.lambda(j)*U1(:,j)*U2(:,j)';
    imagesc(angComp)
    title(['t= ',num2str(timeCoM(j))])
end

figure(3)
U3sort = U3;
for i = 1:R
    subplot(5,4,i)
    j = ind(i);
    plot(U3(:,j))
    U3sort(:,i) = U3(:,j);
end

figure(4)
imagesc((U3sort'))

figure(5)
for i = 1:R
    subplot(5,4,i)
    j = ind(i);
    k = floor(timeCoM(j));
    imagesc(S0(:,:,k))
    title(['t= ',num2str(k)])
end

%% View learned dictionaries
plotDictionary(D1)
plotDictionary(D2)

Y1 = recon(D1,X1);
Y2 = recon(D2,X2);

% Error in time
rel_err = zeros(R,1);
sparsity = zeros(R,1);
dataNorms = zeros(R,1);
vdfs = squeeze(squeeze(sum(sum(X1,1),2)));
for i = 1:R
   rel_err(i) = norm(Y1(:,i)-CP.U{1}(:,i))/norm(CP.U{1}(:,i));
   sparsity(i) = sum(X1(:,:,:,i)>0,'all');
   dataNorms(i) = norm(CP.U{1}(:,i));
   vdfs(:,i) = vdfs(:,i)/max(vdfs(:,i));
end

figure;
subplot(4,1,1)
plot(rel_err)
title('Relative Error')

subplot(4,1,2)
plot(sparsity)
title('Number of Nonzeros')

subplot(4,1,3)
plot(dataNorms)
title('||Data||_2')

subplot(4,1,4)
imagesc(vdfs)
title('VDFs')
xlabel('time')


%% View reconstructions
figure;
subplot(1,2,1)
imagesc(Y1);
set(gca,'visible','off')
title('Data')

subplot(1,2,2)
imagesc(CP.U{1});
set(gca,'visible','off')
title('Recon')

figure;
subplot(1,2,1)
imagesc(Y2);
set(gca,'visible','off')
title('Data')

subplot(1,2,2)
imagesc(CP.U{2});
set(gca,'visible','off')
title('Recon')

function Y = recon(D,X)
Df = fft2(D,1,size(X,2));
Xf = fft2(X);
Yf = sum(bsxfun(@times,Df,Xf),3);
Y = squeeze(real(ifft2(Yf)));
end