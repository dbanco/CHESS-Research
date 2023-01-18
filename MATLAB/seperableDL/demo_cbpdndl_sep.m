
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


% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Construct initial dictionary
D0v = zeros(12,1,144, 'single');
D0v(4:9,:,:) = single(randn(6,1,144));

D0h = zeros(1,12,144, 'single');
D0h(:,4:9,:) = single(randn(1,6,144));
  
for c=1:size(D0h,3)
   D0(:,:,c) = D0v(:,:,c)*D0h(:,:,c);
end

% Set up cbpdndl parameters
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
