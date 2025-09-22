% opt
% Set up algorithm parameters
opt.plotDict = 0;
opt.Verbose = 1;
opt.MaxMainIter = 300;
opt.MaxCGIter = 100;
opt.CGTol = 1e-6;
opt.MaxCGIterX = 100;
opt.CGTolX = 1e-6;
% Rho and sigma params
opt.rho = 1000;
opt.sigma = 500;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.a_min = 1e-4;
opt.lambda_min = 1e-4;
opt.adapt_a = false;
opt.adapt_lambda = false;
opt.NonNegCoef = 1;
opt.NonnegativeDict = 1;
opt.useGpu = false;
opt.Xfixed = 0;
opt.Dfixed = 0;
opt.Recenter = 0;
opt.a = 1;
opt.useMin = false;
opt.AdaptIters = 100;
opt.a_via_lam = true;
opt.l1_iters = 10;
opt.mcdl_init = false;

opt.lambda = 0.1;
opt.lambda2 = 0;
opt.lambda3 = 0;

opt.Penalty = 'l1-norm';
opt.coefInit = 'zeros';
opt.dictInit = 'flat';
opt.Dfixed = 0;
opt.Xfixed = 0;
opt.Recenter = 0;

opt.ZeroMean = 0;

% Regularizer
opt.regularizer = 'filter1';

% Comparison script
rng(1)
[y,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(0.009759000729485,'pseudo-voigt_unmatched');
y = reshape(y,[1,N,T]);

if opt.a_via_lam
    opt.a = 0.95./opt.lambda;
end

% Mean removal and normalisation projections
Pzmn = @(x) bsxfun(@minus, x, mean(mean(x,1),2));
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

% Projection of filter to full image size and its transpose
% (zero-pad and crop respectively)
dsz = [1 M];
Pzp = @(x) x;%zpad(x, xsz(1:2));
PzpT = @(x) bcrop(x, dsz);
Pnonneg = @(x) makeNonneg(x);
 
% Projection of dictionary filters onto constraint set
if opt.ZeroMean
  Pcn = @(x) Pnrm(Pzp(Pzmn(PzpT(x))));
elseif opt.NonnegativeDict
  Pcn = @(x) removeNan(Pnrm(Pzp(PzpT(Pnonneg(x)))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
end

K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);

opt = initXD(opt,N,M,K,J,T,Xtrue,Dtrue);
G = opt.G0;
Gprv = G;
D = Pnrm(opt.G0);
H = zeros(size(G));
X = opt.Y0;
Y = opt.Y0;
U = opt.Y0;
Yprv = Y;

Y = rand(size(Y));
Yf = fft2(Y);

S = reshape(y, [size(y,1) size(y,2) 1 T]);
S = padarray(S,[0 M-1 0 0],0,'pre');
Sfpad = fft2(S);

centerM = round((M+1)/2);
[AG,NormVals,Shifts] = reSampleCustomArrayCenter3(N,G,scales,centerM);
AG = padarray(AG,[0 M-1 0 0],0,'post');
AGf = fft2(AG);
AGSf = bsxfun(@times, conj(AGf), Sfpad);
b = AGSf + opt.rho*fft2(Y-U); % max = 47.902451083396571, 47.902451083396585

% Old
[Xf, ~] = solvemdbi_cg_OF_gpu_zpad(AGf,opt.rho,b,...
            opt.CGTolX, opt.MaxCGIterX, Yf(:),N,M,K,J,T,opt.lambda2,[],[],opt.useGpu,false); 
X1 = ifft2(Xf, 'symmetric');

% New
[Xf, ~] = x_update_switch(Yf(:),AGf,b,opt,N,M,K,J,T);
X2 = ifft2(Xf, 'symmetric');

figure
subplot(2,1,1)
plot(X1(1,:,1,1))
subplot(2,1,2)
plot(X2(1,:,1,1))

norm(X1(:)-X2(:))

% Functions
function u = bcrop(v, sz)
    if numel(sz) <= 2
        if isscalar(sz)
          cs = [sz sz];
        else
          cs = sz;
        end
        u = v(1:cs(1), 1:cs(2), :);
    else
        cs = max(sz,[],2);
        u = zeros(cs(1), cs(2), size(v,3), class(v));
        for k = 1:size(v,3)
            u(1:sz(1,k), 1:sz(2,k), k) = v(1:sz(1,k), 1:sz(2,k), k);
        end
    end
end

function y = makeNonneg(x)
    y=x;
    y(x<0)=0;
end

function y = removeNan(x)
    y=x;
    y(isnan(x))=0;
end
