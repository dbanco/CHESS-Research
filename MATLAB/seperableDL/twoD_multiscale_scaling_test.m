
% Training images
load FlickrCC2_512_512
S0=single(S(:,:,1:20))/255;
nImage = size(S0,3);


%Reduce images size to speed up demo script
tmp = zeros(256, 256, 5, 'single');
for k = 1:size(S0,3),
  tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
end
S0 = tmp;


%% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Construct initial dictionary
N1 = 32;
N2 = 32;
c1 = 2;
c2 = 1;
Ufactors = 3;
K = 1;

D0v = zeros(16,1,K, 'single');
D0v = gaussian_basis_wrap_1D(16,16/2,1.5,'2-norm');
AD0v = reshape(reSampleNu2d(N1,D0v,c1,c2,Ufactors),[N1,1,K,Ufactors]);

D0h = zeros(1,16,K, 'single');
D0h = gaussian_basis_wrap_1D(16,16/2,3,'2-norm');
AD0h = reshape(reSampleNu2d(N2,D0h,c1,c2,Ufactors),[1,N2,K,Ufactors]);

AD0 = zeros(N1,N2,K,Ufactors);

ii= 1;
figure;
for k = 1:K
    for u = 1:Ufactors
        for u2 = 1:Ufactors
        AD0(:,:,k,u) = AD0v(:,:,k,u)*AD0h(:,:,k,u2);
        subplot(Ufactors,Ufactors,ii)
        imagesc(AD0(:,:,k,u))
        set(gca, 'YtickLabel','')
        set(gca, 'XtickLabel','')
        ii = ii + 1;
        end
    end
end

figure;
for k = 1:K
    for u = 1:3
        subplot(1,Ufactors,u)
        plot(AD0v(:,k,u))
        set(gca, 'YtickLabel','')
        set(gca, 'XtickLabel','')
        ii = ii + 1;
    end
end
figure;
for k = 1:K
    for u = 1:3
        subplot(1,Ufactors,u)
        plot(AD0h(k,:,u))
        set(gca, 'YtickLabel','')
        set(gca, 'XtickLabel','')
        ii = ii + 1;
    end
end

figure;
for k = 1:K
    for u = 1:Ufactors
        for u2 = 1:Ufactors
        AD0(:,:,k,u) = AD0v(:,:,k,u)*AD0h(:,:,k,u2);
        subplot(Ufactors,Ufactors,ii)
        imagesc(AD0(:,:,k,u))
        set(gca, 'YtickLabel','')
        set(gca, 'XtickLabel','')
        ii = ii + 1;
        end
    end
end

%

%% Set up cbpdndl parameters
lambda = 0.2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = (50*lambda + 0.5);
opt.sigma = (size(Sh,3));
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.LinSolve = 'CG';
opt.CGTol = 1e-3;
opt.CGTolAuto = 0;

% Do dictionary learning

[D, X, Dx, Dy, optinf] = cbpdndl_sep(D0, D0h, D0v, Sh, lambda, opt);

saveOutputsCBPDNDL_Sep(topDir,ending,Y,D, X, Dx, Dy, optinf,lambda)

% Display learned dictionary
figure;
imdisp(tiledict(D));

% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');
