function [Y,Ymin,Uvel,Vvel,optinf, Jfn, relErr] = cbpdn_cg_OF_multiScales_gpu(D0, S, lambda, lambda2, opt, scales,Uvel,Vvel)
% cbpdndl -- Convolutional BPDN Dictionary Learning
%
%         argmin_{x_m,d_m} (1/2) \sum_k ||\sum_m d_m * x_k,m - s_k||_2^2 +
%                           lambda \sum_k \sum_m ||x_k,m||_1
%
%         The solution is computed using Augmented Lagrangian methods
%         (see boyd-2010-distributed) with efficient solution of the
%         main linear systems (see wohlberg-2016-efficient).
%
% Usage:
%       [D, Y, optinf] = cbpdndl(D0, S, lambda, opt)
%
% Input:
%       D0          Initial dictionary
%       S           Input images
%       lambda      Regularization parameter
%       opt         Options/algorithm parameters structure (see below)
%
% Output:
%       D           Dictionary filter set (3D array)
%       X           Coefficient maps (4D array)
%       optinf      Details of optimisation
%
%
% Options structure fields:
%   Verbose          Flag determining whether iteration status is displayed.
%                    Fields are iteration number, functional value,
%                    data fidelity term, l1 regularisation term, and
%                    primal and dual residuals (see Sec. 3.3 of
%                    boyd-2010-distributed). The values of rho and sigma
%                    are also displayed if options request that they are
%                    automatically adjusted.
%   MaxMainIter      Maximum main iterations
%   AbsStopTol       Absolute convergence tolerance (see Sec. 3.3.1 of
%                    boyd-2010-distributed)
%   RelStopTol       Relative convergence tolerance (see Sec. 3.3.1 of
%                    boyd-2010-distributed)
%   L1Weight         Weight array for L1 norm
%   Y0               Initial value for Y
%   U0               Initial value for U
%   G0               Initial value for G (overrides D0 if specified)
%   H0               Initial value for H
%   rho              Augmented Lagrangian penalty parameter
%   AutoRho          Flag determining whether rho is automatically updated
%                    (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoRhoPeriod    Iteration period on which rho is updated
%   RhoRsdlRatio     Primal/dual residual ratio in rho update test
%   RhoScaling       Multiplier applied to rho when updated
%   AutoRhoScaling   Flag determining whether RhoScaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, RhoScaling specifies a maximum allowed
%                    multiplier instead of a fixed multiplier
%   sigma            Augmented Lagrangian penalty parameter
%   AutoSigma        Flag determining whether sigma is automatically
%                    updated (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoSigmaPeriod  Iteration period on which sigma is updated
%   SigmaRsdlRatio   Primal/dual residual ratio in sigma update test
%   SigmaScaling     Multiplier applied to sigma when updated
%   AutoSigmaScaling Flag determining whether SigmaScaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, SigmaScaling specifies a maximum allowed
%                    multiplier instead of a fixed multiplier.
%   StdResiduals     Flag determining whether standard residual definitions
%                    (see Sec 3.3 of boyd-2010-distributed) are used instead
%                    of normalised residuals (see wohlberg-2015-adaptive)
%   XRelaxParam      Relaxation parameter (see Sec. 3.4.3 of
%                    boyd-2010-distributed) for X update
%   DRelaxParam      Relaxation parameter (see Sec. 3.4.3 of
%                    boyd-2010-distributed) for D update
%   LinSolve         Linear solver for main problem: 'SM' or 'CG'
%   MaxCGIter        Maximum CG iterations when using CG solver
%   CGTol            CG tolerance when using CG solver
%   CGTolAuto        Flag determining use of automatic CG tolerance
%   CGTolFactor      Factor by which primal residual is divided to obtain CG
%                    tolerance, when automatic tolerance is active
%   NoBndryCross     Flag indicating whether all solution coefficients
%                    corresponding to filters crossing the image boundary
%                    should be forced to zero.
%   DictFilterSizes  Array of size 2 x M where each column specifies the
%                    filter size (rows x columns) of the corresponding
%                    dictionary filter
%   NonNegCoef       Flag indicating whether solution should be forced to
%                    be non-negative
%   ZeroMean         Force learned dictionary entries to be zero-mean
%
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-12-18
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'License' file distributed with
% the library.


if nargin < 4,
  opt = [];
end
checkopt(opt, defaultopts([]));
opt = defaultopts(opt);

% Set up status display for verbose operation
hstr = ['Itn   Fnc       DFid      l1        OF        HS        XYU       CGIters      '...
        'r(X)      s(X)      '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %4d %9.2e %9.2e';
nsep = 84;
if opt.AutoRho,
  hstr = [hstr '     rho'];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end

sfms = [sfms,'\n'];
if opt.Verbose && opt.MaxMainIter > 0,
  disp(hstr);
  disp(char('-' * ones(1,nsep)));
end

% Collapsing of trailing singleton dimensions greatly complicates
% handling of both SMV and MMV cases. The simplest approach would be
% if S could always be reshaped to 4d, with dimensions consisting of
% image rows, image cols, a single dimensional placeholder for number
% of filters, and number of measurements, but in the single
% measurement case the third dimension is collapsed so that the array
% is only 3d.
switch numel(size(S))
    case 4
        [N1,N2,~,T] = size(S);
    case 3
        [N1,N2,T] = size(S);
    otherwise
        error('Check size of S')
end
M = size(D0,2);
K = size(D0,3);
KJ = 0;
for i = 1:numel(scales)
    KJ = KJ + size(scales{i},2);
end
J = KJ/K;
xsz = [size(S,1) size(S,2) KJ T];
% Insert singleton 3rd dimension (for number of filters) so that
% 4th dimension is number of images in input s volume
S = reshape(S, [size(S,1) size(S,2) 1 T]);

Nx = prod(xsz);
Nd = prod(xsz(1:2))*size(D0,3);
cgt = opt.CGTol;

% Dictionary size may be specified when learning multiscale
% dictionary
if isempty(opt.DictFilterSizes)
  dsz = [size(D0,1) size(D0,2)];
else
  dsz = opt.DictFilterSizes;
end

% Mean removal and normalisation projections
Pzmn = @(x) bsxfun(@minus, x, mean(mean(x,1),2));
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

% Projection of filter to full image size and its transpose
% (zero-pad and crop respectively)
Pzp = @(x) x;%zpad(x, xsz(1:2));
PzpT = @(x) bcrop(x, dsz);
Pnonneg = @(x) makeNonneg(x);
 
% Projection of dictionary filters onto constraint set
if opt.ZeroMean,
  Pcn = @(x) Pnrm(Pzp(Pzmn(PzpT(x))));
elseif opt.NonnegativeDict,
  Pcn = @(x) removeNan(Pnrm(Pzp(PzpT(Pnonneg(x)))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
end

% Start timer
tstart = tic;

% Project initial dictionary onto constraint set
D = Pnrm(D0);

% Compute signal in DFT domain
Sf = fft2(S);

% Set up algorithm parameters and initialise variables
rho = opt.rho;
if isempty(rho), rho = 50*lambda+1; end;
if opt.AutoRho,
  asgr = opt.RhoRsdlRatio;
  asgm = opt.RhoScaling;
end

optinf = struct('itstat', [], 'opt', opt);
rx = Inf;
sx = Inf;
eprix = 0;
eduax = 0;


% Initialise main working variables
if isempty(opt.Y0),
  Y = zeros(xsz, class(S));
else
  Y = opt.Y0;
end
Yprv = Y;
if isempty(opt.U0),
  if isempty(opt.Y0),
    U = zeros(xsz, class(S));
  else
    U = (lambda/rho)*sign(Y);
  end
else
  U = opt.U0;
end


[AG,NormVals] = reSampleCustomArray(N2,D,scales);
AGf = fft2(AG);
AGSf = bsxfun(@times, conj(AGf), Sf);

Yf = fft2(Y);

% Initial solution
k=0;
tk = toc(tstart);

recon = sum(bsxfun(@times,AGf,Yf),3);
Jdf = sum(vec(abs(recon-Sf).^2));
Jl1 = sum(abs(vec(bsxfun(@times, opt.L1Weight, Y))));
Jlg1 = rho*sum((U(:)).^2);
if lambda2 > 0
    [~,~,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(Y,K,Uvel,Vvel,opt.Smoothness,opt.HSiters);
    [Jof, Jhs] = HSobjectivePaper(Fx,Fy,Ft,Uvel,Vvel,K,opt.Smoothness);
else
    Jhs = 0;
    Jof = 0;
end
Jfn = Jdf + lambda*Jl1 + lambda2*Jof + opt.Smoothness*Jhs + Jlg1;

% Initial min solution
Ymin = Y;
minJfn = Jfn;
minCount = 0;

optinf.itstat = [optinf.itstat;...
       [k Jfn Jdf Jl1 Jof Jhs Jlg1 rx sx eprix eduax rho tk]];
  if opt.Verbose,
    dvc = [k, Jfn, Jdf, Jl1, Jof, Jhs, Jlg1, 0,0, rx, sx];
    if opt.AutoRho,
      dvc = [dvc rho];
    end
    fprintf(sfms, dvc);
  end

% Main loop
k = 1;
while k <= opt.MaxMainIter && (rx > eprix|sx > eduax),
    % Solve X subproblem. It would be simpler and more efficient (since the
    % DFT is already available) to solve for X using the main dictionary
    % variable D as the dictionary, but this appears to be unstable. Instead,
    % use the projected dictionary variable G

    % Solve subproblem
%     recon = sum(bsxfun(@times,AGf,Yf),3);
%     Jdf1 = sum(vec(abs(recon-Sf).^2))
    [Xf, cgst] = solvemdbi_cg_OF_gpu(AGf, rho, AGSf + rho*fft2(Y - U) ,...
        opt.CGTolX, opt.MaxCGIterX, Yf(:),N2,K,J,T,lambda2,Uvel,Vvel); 
    cgIters2 = cgst.pit;
    X = ifft2(Xf, 'symmetric');

    % Data fidelity term in Fourier domain
%     recon = sum(bsxfun(@times,AGf,Xf),3);
%     Jdf2 = sum(vec(abs(recon-Sf).^2))
    


    clear Xf;
    
    % See pg. 21 of boyd-2010-distributed
    if opt.XRelaxParam == 1,
        Xr = X;
    else
        Xr = opt.XRelaxParam*X + (1-opt.XRelaxParam)*Y;
    end
    
    % Solve Y subproblem
    Y = shrink(Xr + U, (lambda/rho)*opt.L1Weight);
    if opt.NonNegCoef,
        Y(Y < 0) = 0;
    end
    
    if opt.NoBndryCross,
        %Y((end-max(dsz(1,:))+2):end,:,:,:) = 0;
        Y((end-size(D0,1)+2):end,:,:,:) = 0;
        %Y(:,(end-max(dsz(2,:))+2):end,:,:) = 0;
        Y(:,(end-size(D0,2)+2):end,:,:) = 0;
    end
    Yf = fft2(Y);
    
    % Data fidelity term in Fourier domain
%     recon = sum(bsxfun(@times,AGf,Yf),3);
%     Jdf3 = sum(vec(abs(recon-Sf).^2))

    % Update dual variable corresponding to X, Y, Z
    U = U + Xr - Y;
    clear DPXr Xr;
    
    % Compute primal and dual residuals and stopping thresholds for X update
    nX = norm(X(:)); nY = norm(Y(:)); nU = norm(U(:)); 
    if opt.StdResiduals,
        % See pp. 19-20 of boyd-2010-distributed
        rx = norm(vec(X - Y));
        sx = norm(vec(rho*(Yprv - Y)));
        eprix = sqrt(Nx)*opt.AbsStopTol+max(nX,nY)*opt.RelStopTol;
        eduax = sqrt(Nx)*opt.AbsStopTol+rho*nU*opt.RelStopTol;
    else
        % See wohlberg-2015-adaptive
        rx = norm(vec(X - Y))/max(nX,nY);
        sx = norm(vec(Yprv - Y))/nU;
        eprix = sqrt(Nx)*opt.AbsStopTol/max(nX,nY)+opt.RelStopTol;
        eduax = sqrt(Nx)*opt.AbsStopTol/(rho*nU)+opt.RelStopTol;
    end

    Jlg1 = rho*sum((X(:)-Y(:)+U(:)).^2);

    clear X;
    
    % Update record of previous step Y
    Yprv = Y;
    
    % Apply CG auto tolerance policy if enabled
    if opt.CGTolAuto && (rd/opt.CGTolFactor) < cgt,
        cgt = rd/opt.CGTolFactor;
    end

    % Update optical flow velocities
    if (opt.UpdateVelocity && (lambda2 > 0)) || nargin < 6
        [Uvel,Vvel,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(Y,K,Uvel,Vvel,opt.Smoothness,opt.HSiters);
    end

    % Data fidelity term in Fourier domain
    recon = sum(bsxfun(@times,AGf,Yf),3);
    Jdf = sum(vec(abs(recon-Sf).^2));
    
    % Sparsity term
    Jl1 = sum(abs(vec(bsxfun(@times, opt.L1Weight, Y))));
    
    % Optical flow terms
    if lambda2 > 0
        [Jof, Jhs] = HSobjectivePaper(Fx,Fy,Ft,Uvel,Vvel,K,opt.Smoothness);
    else
        Jof = 0;
        Jhs = 0;
    end
    
    % Full objective
    Jfn = Jdf + lambda*Jl1 + lambda2*Jof + opt.Smoothness*Jhs + Jlg1;

    if Jfn < minJfn
        Ymin = Y;
        minJfn = Jfn;
        minCount = 0;
    else
        minCount = minCount + 1;
        if minCount >10
        %           break
        end
    end
    
    % Record and display iteration details
    tk = toc(tstart);
    optinf.itstat = [optinf.itstat;...
    [k Jfn Jdf Jl1 Jof Jhs Jlg1 rx sx eprix eduax rho tk]];
    if opt.Verbose
        dvc = [k, Jfn, Jdf, Jl1, Jof, Jhs, Jlg1,cgIters2, rx, sx];
        if opt.AutoRho
            dvc = [dvc rho];
        end
        fprintf(sfms, dvc);
    end

  % See wohlberg-2015-adaptive and pp. 20-21 of boyd-2010-distributed
  if opt.AutoRho
    if k ~= 1 && mod(k, opt.AutoRhoPeriod) == 0
      if opt.AutoRhoScaling
        rhomlt = sqrt(rx/sx);
        if rhomlt < 1, rhomlt = 1/rhomlt; end
        if rhomlt > opt.RhoScaling, rhomlt = opt.RhoScaling; end
      else
        rhomlt = opt.RhoScaling;
      end
      rsf = 1;
      if rx > opt.RhoRsdlRatio*sx, rsf = rhomlt; end
      if sx > opt.RhoRsdlRatio*rx, rsf = 1/rhomlt; end
      rho = rsf*rho;
      U = U/rsf;
    end
  end
  
  k = k + 1;

end

relErr = Jdf/sum(vec((Sf).^2));

% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Y = Y;
optinf.U = U;

optinf.lambda = lambda;
optinf.rho = rho;
optinf.cgt = cgt;
if exist('cgst'), optinf.cgst = cgst; end

if opt.Verbose && opt.MaxMainIter > 0,
  disp(char('-' * ones(1,nsep)));
end

if opt.plotDict
    close(fdict);
end

return


function u = vec(v)

  u = v(:);

return


function u = shrink(v, lambda)

  if isscalar(lambda),
    u = sign(v).*max(0, abs(v) - lambda);
  else
    u = sign(v).*max(0, bsxfun(@minus, abs(v), lambda));
  end

return


function u = zpad(v, sz)

  u = zeros(sz(1), sz(2), size(v,3), size(v,4), class(v));
  u(1:size(v,1), 1:size(v,2),:,:) = v;

return


function u = bcrop(v, sz)

  if numel(sz) <= 2,
    if numel(sz) == 1
      cs = [sz sz];
    else
      cs = sz;
    end
    u = v(1:cs(1), 1:cs(2), :);
  else
    cs = max(sz,[],2);
    u = zeros(cs(1), cs(2), size(v,3), class(v));
    for k = 1:size(v,3),
      u(1:sz(1,k), 1:sz(2,k), k) = v(1:sz(1,k), 1:sz(2,k), k);
    end
  end

return

function y = makeNonneg(x)
    y=x;
    y(x<0)=0;
return

function y = removeNan(x)
    y=x;
    y(isnan(x))=0;
return

function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'AbsStopTol'),
    opt.AbsStopTol = 1e-5;
  end
  if ~isfield(opt,'RelStopTol'),
    opt.RelStopTol = 1e-3;
  end
  if ~isfield(opt,'L1Weight'),
    opt.L1Weight = 1;
  end
  if ~isfield(opt,'Y0'),
    opt.Y0 = [];
  end
  if ~isfield(opt,'U0'),
    opt.U0 = [];
  end
  if ~isfield(opt,'Z0'),
    opt.Z0 = [];
  end
  if ~isfield(opt,'V0'),
    opt.V0 = [];
  end
  if ~isfield(opt,'G0'),
    opt.G0 = [];
  end
  if ~isfield(opt,'H0'),
    opt.H0 = [];
  end
  if ~isfield(opt,'rho'),
    opt.rho = [];
  end
  if ~isfield(opt,'AutoRho'),
    opt.AutoRho = 0;
  end
  if ~isfield(opt,'AutoRhoPeriod'),
    opt.AutoRhoPeriod = 10;
  end
  if ~isfield(opt,'RhoRsdlRatio'),
    opt.RhoRsdlRatio = 10;
  end
  if ~isfield(opt,'RhoScaling'),
    opt.RhoScaling = 2;
  end
  if ~isfield(opt,'AutoRhoScaling'),
    opt.AutoRhoScaling = 0;
  end
  if ~isfield(opt,'sigma'),
    opt.sigma = [];
  end
  if ~isfield(opt,'AutoSigma'),
    opt.AutoSigma = 0;
  end
  if ~isfield(opt,'AutoSigmaPeriod'),
    opt.AutoSigmaPeriod = 10;
  end
  if ~isfield(opt,'SigmaRsdlRatio'),
    opt.SigmaRsdlRatio = 10;
  end
  if ~isfield(opt,'SigmaScaling'),
    opt.SigmaScaling = 2;
  end
  if ~isfield(opt,'AutoSigmaScaling'),
    opt.AutoSigmaScaling = 0;
  end
  if ~isfield(opt,'StdResiduals'),
    opt.StdResiduals = 0;
  end
  if ~isfield(opt,'XRelaxParam'),
    opt.XRelaxParam = 1;
  end
  if ~isfield(opt,'DRelaxParam'),
    opt.DRelaxParam = 1;
  end
  if ~isfield(opt,'LinSolve'),
    opt.LinSolve = 'SM';
  end
  if ~isfield(opt,'MaxCGIter'),
    opt.MaxCGIter = 1000;
  end
  if ~isfield(opt,'CGTol'),
    opt.CGTol = 1e-3;
  end
  if ~isfield(opt,'CGTolX'),
    opt.CGTolX = 1e-3;
  end
  if ~isfield(opt,'CGTolAuto'),
    opt.CGTolAuto = 0;
  end
  if ~isfield(opt,'CGTolAutoFactor'),
    opt.CGTolFactor = 50;
  end
  if ~isfield(opt,'NoBndryCross'),
    opt.NoBndryCross = 0;
  end
  if ~isfield(opt,'DictFilterSizes'),
    opt.DictFilterSizes = [];
  end
  if ~isfield(opt,'NonNegCoef'),
    opt.NonNegCoef = 0;
  end
  if ~isfield(opt,'ZeroMean'),
    opt.ZeroMean = 0;
  end
  if ~isfield(opt,'NonnegativeDict'),
    opt.NonnegativeDict = 0;
  end
  if ~isfield(opt,'plotDict'),
    opt.plotDict = 0;
  end
return
