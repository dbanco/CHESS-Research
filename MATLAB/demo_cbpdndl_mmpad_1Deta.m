% Script demonstrating usage of the cbpdndl function.
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-07-30
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.

% Training images
% Load MMPAD sequence
r = 1;
data_dir = ['D:\MMPAD_data_nr1\ring', num2str(r), '_zero'];
T = 200;

for t = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t),'.mat']))
    if t == 1
        [N,M] = size(polar_image);
        S0 = zeros(1,M,T);
    end
    S0(1,:,t) = sum(polar_image,1);
end
S0 = S0(1,:,1:4:T)/1e3;

% Set up cbpdndl parameters
lambda = 50;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 250;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(S0,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.DictFilterSizes = [ones(1,15);
                       12,12,12,12,12,20,20,20,20,40,40,40,64,64,88];
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;

P.var_theta = logspace(-0.3,1.699,15).^2;
P.num_theta = 88;
P.mean = 6;
P.basis = 'norm2';
D0 = zeros(1,88,15);
for i = 1:15
    D0(1,:,:) = dictionary(P);
end


% Do dictionary learning
[D, X, optinf] = cbpdndl(D0, S0, lambda, opt);

% Display learned dictionary
figure;
imagesc(squeeze(D));

% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');

%% Load reconstructions
% load('dict_learn_results_1D.mat')
Df = fft2(D,size(S0,1),size(S0,2));
% X(X<0.1) = 0;
Xf = fft2(X);
Yf = sum(bsxfun(@times,Df,Xf),3);
Y = squeeze(real(ifft2(Yf)));

%% Display
% tt = 37;
% Arrange images
dispY = arrangeIms(Y);
dispS = arrangeIms(S0);
dispD = arrangeIms(D);

% figure;
% imagesc(squeeze(D));
plotDictionary(D)



figure;
subplot(1,2,1)
imagesc(dispY);
set(gca,'visible','off')
title('Data')

subplot(1,2,2)
imagesc(dispS');
set(gca,'visible','off')
title('Recon')

% Error in time
rel_err = zeros(50,1);
sparsity = zeros(50,1);
dataNorms = zeros(50,1);
vdfs = squeeze(squeeze(sum(sum(X,1),2)));
for t = 1:50
   rel_err(t) = norm(Y(:,t)'-S0(1,:,t))/norm(S0(1,:,t));
   sparsity(t) = sum(X(:,:,:,t)>0,'all');
   dataNorms(t) = norm(S0(1,:,t));
   vdfs(:,t) = vdfs(:,t)/max(vdfs(:,t));
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
