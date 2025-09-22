function [D, Y, X, Dmin, Ymin,optinf, Jfn, relErr] = mcdl_regularized_gpu(D0, S, opt, scales)
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
%                    boyd-2010-distributed). The values of rho1 and rho2
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
%   rho1              Augmented Lagrangian penalty parameter
%   AutoRho1          Flag determining whether rho1 is automatically updated
%                    (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoRho1Period    Iteration period on which rho1 is updated
%   Rho1RsdlRatio     Primal/dual residual ratio in rho1 update test
%   Rho1Scaling       Multiplier applied to rho1 when updated
%   AutoRho1Scaling   Flag determining whether Rho1Scaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, Rho1Scaling specifies a maximum allowed
%                    multiplier instead of a fixed multiplier
%   rho2            Augmented Lagrangian penalty parameter
%   AutoRho2        Flag determining whether rho2 is automatically
%                    updated (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoRho2Period  Iteration period on which rho2 is updated
%   Rho2RsdlRatio   Primal/dual residual ratio in rho2 update test
%   Rho2Scaling     Multiplier applied to rho2 when updated
%   AutoRho2Scaling Flag determining whether Rho2Scaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, Rho2Scaling specifies a maximum allowed
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
%   Recenter         Recenters learned dictionary atom before sparse coding
%
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-12-18
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'License' file distributed with
% the library.


if nargin < 4
  opt = [];
end
checkopt(opt, defaultopts([]));
opt = defaultopts(opt);

% Set up status display for verbose operation
hstr = ['Itn   Fnc       DFid      l1        Reg      Jrho1     Jrho2     Cnstr     CGIters      '...
        'r(X)      s(X)      r(D)      s(D) '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %4d %4d %9.2e %9.2e %9.2e %9.2e';
nsep = 84;
if opt.AutoRho1
  hstr = [hstr '     rho1'];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
if opt.AutoRho2
  hstr = [hstr '     rho2  '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
sfms = [sfms,'\n'];
if opt.Verbose && opt.MaxMainIter > 0
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
        [~,N2,~,T] = size(S);
    case 3
        [~,N2,T] = size(S);
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
centerM = round((M+1)/2);
xsz = [size(S,1) size(S,2) KJ T];
% Insert singleton 3rd dimension (for number of filters) so that
% 4th dimension is number of images in input s volume
S = reshape(S, [size(S,1) size(S,2) 1 T]);
Spad = padarray(S,[0 M-1 0 0],0,'pre');

Nx = prod(xsz);
Nd = prod(xsz(1:2))*size(D0,3);
cgt = opt.CGTol;

% Dictionary size
if isempty(opt.DictFilterSizes)
  dsz = [size(D0,1) size(D0,2)];
else
  dsz = opt.DictFilterSizes;
end

lambda = opt.lambda;
lambda2 = opt.lambda2;
% Set a paramter via lambda
if opt.a_via_lam
    opt.a = 0.95./lambda;
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
if opt.ZeroMean
  Pcn = @(x) Pnrm(Pzp(Pzmn(PzpT(x))));
elseif opt.NonnegativeDict
  Pcn = @(x) removeNan(Pnrm(Pzp(PzpT(Pnonneg(x)))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
end

% Start timer
tstart = tic;

% Project initial dictionary onto constraint set
D = Pnrm(D0);


% Compute signal in DFT domain
Sfpad = fft2(Spad);

% Set up algorithm parameters and initialise variables
rho1 = opt.rho1;
if isempty(rho1), rho1 = 50*lambda+1; end
% if opt.AutoRho1
%   asgr = opt.Rho1RsdlRatio;
%   asgm = opt.Rho1Scaling;
% end
rho2 = opt.rho2;
if isempty(rho2), rho2 = size(S,3); end
% if opt.AutoRho2
%   asdr = opt.Rho2RsdlRatio;
%   asdm = opt.Rho2Scaling;
% end
optinf = struct('itstat', [], 'opt', opt);
rx = Inf;
sx = Inf;
rd = Inf;
sd = Inf;
eprix = 0;
eduax = 0;
eprid = 0;
eduad = 0;

% Initialise main working variables
if isempty(opt.Y0)
  Y = zeros(xsz, class(S));
else
  Y = opt.Y0;
end
Yprv = Y;
X = Y;

if isempty(opt.U0)
  if isempty(opt.Y0)
    U = zeros(xsz, class(S));
  else
    U = (lambda/rho1)*sign(Y);
  end
else
  U = opt.U0;
end

if isempty(opt.G0)
  G = Pzp(D);
else
  G = opt.G0;
end
Gprv = G;
if isempty(opt.H0)
    H = zeros(size(G), class(S));
else
  H = opt.H0;
end

Dr = D;
Xr = X;

if opt.plotDict
    fdict = figure;
    % frecon = figure;
    freconA = figure;
    xFig = figure;
end

[AG,NormVals,Shifts] = reSampleCustomArrayCenter3(N2,G,scales,centerM);
AGpad = padarray(AG,[0 M-1 0 0],0,'post');

AGf = fft2(AGpad);
AGSf = bsxfun(@times, conj(AGf), Sfpad);
Yf = fft2(Y);

if opt.useGpu
    Atop = @(r) ifft2(pagefun(@times, conj(AGf), r),'symmetric');
else
    Atop = @(r) ifft2(bsxfun(@times, conj(AGf), r),'symmetric');
end

% Initial solution
k=0;
tk = toc(tstart);
Jcn = norm(vec(Pcn(D) - D));
recon = unpad(ifft2(sum(bsxfun(@times,AGf,Yf),3),'symmetric'),M-1,'pre');
Jdf = 0.5*norm(recon-S,'fro').^2;

a_n = opt.a;
Jl1 = lambda*compute_sparse_penalty(Y,a_n,opt,k);
Jrho1 = 0.5*rho1*norm(U,'fro').^2;
Jrho2 = 0.5*rho2*norm(H,'fro').^2;

if lambda2 > 0
    Jreg = lambda2*compute_penalty(Y,K,J,opt.regularizer);
else
    Jreg = 0;
end

Jfn = Jdf + Jl1 + Jreg + Jrho1 + Jrho2 ;

% Initial min solution
Ymin = Y;
Gmin = G;
minJfn = Jfn;
minCount = 0;

optinf.itstat = [optinf.itstat;...
       [k Jfn Jdf Jl1 Jreg Jrho1 Jrho2 rx sx rd sd eprix eduax eprid eduad rho1 rho2 tk]];
if opt.Verbose
    dvc = [k, Jfn, Jdf, Jl1, Jreg, Jrho1, Jrho2, Jcn, 0,0, rx, sx, rd, sd];
    if opt.AutoRho1
        dvc = [dvc rho1];
    end
    if opt.AutoRho2
        dvc = [dvc rho2];
    end
    fprintf(sfms, dvc);
end

if opt.plotDict
    figure(freconA)
    plot(squeeze(S(:,:,:,10)))
    hold on
    plot(squeeze(recon(:,:,:,10)))
    hold off
    
    figure(xFig)
    ii = 1;
    for kk = 1:K
        subplot(1,K,kk)
        Ui = size(scales{kk},2) + ii - 1;
        imagesc( squeeze(sum(Y(:,:,ii:Ui,:),3)) )
        ii = ii + Ui;
    end
    pause()
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Main loop %%%%%%%%%%
k = 1;
while k <= opt.MaxMainIter && (rx > eprix||sx > eduax||rd > eprid||sd >eduad)
    % Solve D subproblem. Similarly, it would be simpler and more efficient to
    % solve for D using the main coefficient variable X as the coefficients,
    % but it appears to be more stable to use the shrunk coefficient variable Y
    
    if ~opt.Dfixed && k > 1
        AYS = reSampleTransCustomArrayCenter3(M,ifft2(sum(bsxfun(@times, conj(Yf), Sfpad), 4),'symmetric'),scales,centerM,NormVals,Shifts);
        [D, cgst] = solvemdbi_cg_multirate_custom_gpu_zpad_center3(Yf, rho2, AYS + rho2*(G - H),...
                          cgt, opt.MaxCGIter, G(:),N2,M,scales,NormVals,Shifts,centerM,opt.useGpu);
        cgIters1 = cgst.pit;
        
        % See pg. 21 of boyd-2010-distributed
        if opt.DRelaxParam == 1
            Dr = D;
        else
            Dr = opt.DRelaxParam*D + (1-opt.DRelaxParam)*G;
        end
        
        % Solve G subproblem
        G = Pcn(Dr + H);

        % Update alphas
        [AG,NormVals,Shifts] = reSampleCustomArrayCenter3(N2,G,scales,centerM);
        AG = padarray(AG,[0 M-1 0 0],0,'post');
        AGf = fft2(AG);
        AGSf = bsxfun(@times, conj(AGf), Sfpad);
        
        if opt.useGpu
            Atop = @(r) ifft2(pagefun(@times, conj(AGf), r),'symmetric');
        else
            Atop = @(r) ifft2(bsxfun(@times, conj(AGf), r),'symmetric');
        end
        
        % Update dual variable corresponding to D, G
        H = H + Dr - G;
    else
        cgIters1 = 0;
    end

    % Compute primal and dual residuals and stopping thresholds for D update
    nD = norm(D(:)); nG = norm(G(:)); nH = norm(H(:));
    if opt.StdResiduals
        % See pp. 19-20 of boyd-2010-distributed
        rd = norm(vec(D - G));
        sd = norm(vec(rho2*(Gprv - G)));
        eprid = sqrt(Nd)*opt.AbsStopTol+max(nD,nG)*opt.RelStopTol;
        eduad = sqrt(Nd)*opt.AbsStopTol+rho2*nH*opt.RelStopTol;
    else
        % See wohlberg-2015-adaptive
        rd = norm(vec(D - G))/max(nD,nG);
        sd = norm(vec(Gprv - G))/nH;
        eprid = sqrt(Nd)*opt.AbsStopTol/max(nD,nG)+opt.RelStopTol;
        eduad = sqrt(Nd)*opt.AbsStopTol/(rho2*nH)+opt.RelStopTol;
    end

    % Solve X subproblem. It would be simpler and more efficient (since the
    % DFT is already available) to solve for X using the main dictionary
    % variable D as the dictionary, but this appears to be unstable. Instead,
    % use the projected dictionary variable G
    if ~opt.Xfixed
        bf = AGSf + rho1*fft2(Y-U);
        [Xf, cgst, opt] = x_update_switch(Y,Yf,AGf,bf,Spad,Y,U,opt,N2,M,K,J,T,Atop);

        cgIters2 = cgst.pit;

        X = ifft2(Xf, 'symmetric');

        clear Xf;
        
        % See pg. 21 of boyd-2010-distributed
        if opt.XRelaxParam == 1
            Xr = X;
        else
            Xr = opt.XRelaxParam*X + (1-opt.XRelaxParam)*Y;
        end
        
        if opt.adapt_a && k <= (opt.AdaptIters + opt.l1_iters)
            a_n = opt.a_min*(opt.a/opt.a_min)^((k-opt.l1_iters)/opt.AdaptIters);
        else
            a_n = opt.a;
        end

        % Solve Y subproblem
        Y = apply_sparse_penalty(Xr + U,lambda,k,opt,a_n);
        Yf = fft2(Y);

        % Update dual variable corresponding to X, Y, Z
        U = U + Xr - Y;
    else
        X = Y;
        cgIters2 = 0;
    end

    % Compute primal and dual residuals and stopping thresholds for X update
    nX = norm(X(:)); nY = norm(Y(:)); nU = norm(U(:)); 
    if opt.StdResiduals
        % See pp. 19-20 of boyd-2010-distributed
        rx = norm(vec(X - Y));
        sx = norm(vec(rho1*(Yprv - Y)));
        eprix = sqrt(Nx)*opt.AbsStopTol+max(nX,nY)*opt.RelStopTol;
        eduax = sqrt(Nx)*opt.AbsStopTol+rho1*nU*opt.RelStopTol;
    else
        % See wohlberg-2015-adaptive
        rx = norm(vec(X - Y))/max(nX,nY);
        sx = norm(vec(Yprv - Y))/nU;
        eprix = sqrt(Nx)*opt.AbsStopTol/max(nX,nY)+opt.RelStopTol;
        eduax = sqrt(Nx)*opt.AbsStopTol/(rho1*nU)+opt.RelStopTol;
    end

    % Update record of previous step Y
    Yprv = Y;
    
    % Apply CG auto tolerance policy if enabled
    if opt.CGTolAuto && (rd/opt.CGTolFactor) < cgt
        cgt = rd/opt.CGTolFactor;
    end
    
    % Compute measure of D constraint violation
    Jcn = norm(vec(Pcn(D) - D));
    
    % Update record of previous step G
    Gprv = G;
    
    % Data fidelity term 
    recon = unpad(ifft2(sum(bsxfun(@times,AGf,Yf),3),'symmetric'),M-1,'pre');
    Jdf = 0.5*norm(recon-S,'fro').^2;

    % Sparsity term
    Jl1 = lambda*compute_sparse_penalty(Y,a_n,opt,k);

    % Regularization term
    if lambda2 > 0
        Jreg = lambda2*compute_penalty(Y,K,J,opt.regularizer);
    else
        Jreg = 0;
    end

    % Constraint terms
    Jrho1 = 0.5*rho1*norm(Xr-Y,'fro').^2;
    Jrho2 = 0.5*rho2*norm(Dr-G,'fro').^2;
    
    % Full objective
    Jfn = Jdf + Jl1 + Jreg + Jrho1 + Jrho2;

    % Plot dictionary progress
    if opt.plotDict
        figure(freconA)
        plot(squeeze(S(:,:,:,10)))
        hold on
        plot(squeeze(recon(:,:,:,10)))
        hold off
        
        figure(xFig)
        ii = 1;
        for kk = 1:K
            subplot(1,K,kk)
            Ui = size(scales{kk},2) + ii - 1;
            imagesc( squeeze(sum(Y(:,:,ii:Ui,:),3)) )
            ii = ii + Ui;
        end
        pause()
    end 

    if opt.plotDict
        AG = reSampleCustomArrayCenter3(N2,G,scales,centerM);
        plotDictUsage(AG,K,1,fdict);
        pause(0.0001)
    end    
    
    if Jfn < minJfn
        Ymin = Y;
        Gmin = G;
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
    [k Jfn Jdf Jl1 Jreg Jrho1 Jrho2 rx sx rd sd eprix eduax eprid eduad rho1 rho2 tk]];
    if opt.Verbose
        if opt.AutoRho1 && opt.AutoRho2
            dvc = [k, Jfn, Jdf, Jl1, Jreg, Jrho1, Jrho2, Jcn,...
                   cgIters1,cgIters2, rx, sx, rd, sd, rho1, rho2];
        elseif opt.AutoRho2
            dvc = [k, Jfn, Jdf, Jl1, Jreg, Jrho1, Jrho2, Jcn,...
                   cgIters1,cgIters2, rx, sx, rd, sd, rho2];
        elseif opt.AutoRho1
            dvc = [k, Jfn, Jdf, Jl1, Jreg, Jrho1, Jrho2, Jcn,...
                   cgIters1,cgIters2, rx, sx, rd, sd, rho1];
        else
            dvc = [k, Jfn, Jdf, Jl1, Jreg, Jrho1, Jrho2, Jcn,...
                   cgIters1,cgIters2, rx, sx, rd, sd];
        end
        fprintf(sfms, dvc);
    end

  % See wohlberg-2015-adaptive and pp. 20-21 of boyd-2010-distributed
  if opt.AutoRho1
    if k ~= 1 && mod(k, opt.AutoRho1Period) == 0
      if opt.AutoRho1Scaling
        rhomlt = sqrt(rx/sx);
        if rhomlt < 1, rhomlt = 1/rhomlt; end
        if rhomlt > opt.Rho1Scaling, rhomlt = opt.Rho1Scaling; end
      else
        rhomlt = opt.Rho1Scaling;
      end
      rsf = 1;
      if rx > opt.Rho1RsdlRatio*sx, rsf = rhomlt; end
      if sx > opt.Rho1RsdlRatio*rx, rsf = 1/rhomlt; end
      rho1 = rsf*rho1;
      opt.rho1 = rho1;
      U = U/rsf;
    end
  end
  if opt.AutoRho2
    if k ~= 1 && mod(k, opt.AutoRho2Period) == 0
      if opt.AutoRho2Scaling
        sigmlt = sqrt(rd/sd);
        if sigmlt < 1, sigmlt = 1/sigmlt; end
        if sigmlt > opt.Rho2Scaling, sigmlt = opt.Rho2Scaling; end
      else
        sigmlt = opt.Rho2Scaling;
      end
      ssf = 1;
      if rd > opt.Rho2RsdlRatio*sd, ssf = sigmlt; end
      if sd > opt.Rho2RsdlRatio*rd, ssf = 1/sigmlt; end
      rho2 = ssf*rho2;
      opt.rho2 = rho2;
      H = H/ssf;
    end
  end

  k = k + 1;

end

D = PzpT(G);
Dmin = PzpT(Gmin);
relErr = Jdf/sum(vec((Sfpad).^2));

% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Y = Y;
optinf.U = U;
optinf.G = G;
optinf.H = H;
optinf.lambda = lambda;
optinf.rho1 = rho1;
optinf.rho2 = rho2;
optinf.cgt = cgt;
if exist('cgst','var'), optinf.cgst = cgst; end

if opt.Verbose && opt.MaxMainIter > 0
  disp(char('-' * ones(1,nsep)));
end

if opt.plotDict
    close(fdict);
end

return

function u = vec(v)

  u = v(:);

return

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
  if ~isfield(opt,'Verbose')
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter')
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'AbsStopTol')
    opt.AbsStopTol = 1e-5;
  end
  if ~isfield(opt,'RelStopTol')
    opt.RelStopTol = 1e-3;
  end
  if ~isfield(opt,'L1Weight')
    opt.L1Weight = 1;
  end
  if ~isfield(opt,'Y0')
    opt.Y0 = [];
  end
  if ~isfield(opt,'U0')
    opt.U0 = [];
  end
  if ~isfield(opt,'Z0')
    opt.Z0 = [];
  end
  if ~isfield(opt,'V0')
    opt.V0 = [];
  end
  if ~isfield(opt,'G0')
    opt.G0 = [];
  end
  if ~isfield(opt,'H0')
    opt.H0 = [];
  end
  if ~isfield(opt,'rho1')
    opt.rho1 = [];
  end
  if ~isfield(opt,'AutoRho1')
    opt.AutoRho1 = 0;
  end
  if ~isfield(opt,'AutoRho1Period')
    opt.AutoRho1Period = 10;
  end
  if ~isfield(opt,'Rho1RsdlRatio')
    opt.Rho1RsdlRatio = 10;
  end
  if ~isfield(opt,'Rho1Scaling')
    opt.Rho1Scaling = 2;
  end
  if ~isfield(opt,'AutoRho1Scaling')
    opt.AutoRho1Scaling = 0;
  end
  if ~isfield(opt,'rho2')
    opt.rho2 = [];
  end
  if ~isfield(opt,'AutoRho2')
    opt.AutoRho2 = 0;
  end
  if ~isfield(opt,'AutoRho2Period')
    opt.AutoRho2Period = 10;
  end
  if ~isfield(opt,'Rho2RsdlRatio')
    opt.Rho2RsdlRatio = 10;
  end
  if ~isfield(opt,'Rho2Scaling')
    opt.Rho2Scaling = 2;
  end
  if ~isfield(opt,'AutoRho2Scaling')
    opt.AutoRho2Scaling = 0;
  end
  if ~isfield(opt,'StdResiduals')
    opt.StdResiduals = 1;
  end
  if ~isfield(opt,'XRelaxParam')
    opt.XRelaxParam = 1;
  end
  if ~isfield(opt,'DRelaxParam')
    opt.DRelaxParam = 1;
  end
  if ~isfield(opt,'LinSolve')
    opt.LinSolve = 'SM';
  end
  if ~isfield(opt,'MaxCGIter')
    opt.MaxCGIter = 1000;
  end
  if ~isfield(opt,'CGTol')
    opt.CGTol = 1e-3;
  end
  if ~isfield(opt,'CGTolX')
    opt.CGTolX = 1e-3;
  end
  if ~isfield(opt,'CGTolAuto')
    opt.CGTolAuto = 0;
  end
  if ~isfield(opt,'CGTolAutoFactor')
    opt.CGTolFactor = 50;
  end
  if ~isfield(opt,'NoBndryCross')
    opt.NoBndryCross = 0;
  end
  if ~isfield(opt,'DictFilterSizes')
    opt.DictFilterSizes = [];
  end
  if ~isfield(opt,'NonNegCoef')
    opt.NonNegCoef = 0;
  end
  if ~isfield(opt,'ZeroMean')
    opt.ZeroMean = 0;
  end
  if ~isfield(opt,'NonnegativeDict')
    opt.NonnegativeDict = 0;
  end
  if ~isfield(opt,'plotDict')
    opt.plotDict = 0;
  end
  if ~isfield(opt,'useGpu')
      opt.useGpu = 1;
  end
  if ~isfield(opt,'Xfixed')
      opt.Xfixed = 0;
  end
  if ~isfield(opt,'Dfixed')
      opt.Dfixed = 0;
  end
  if ~isfield(opt,'Penalty')
      opt.Penalty = 'l1-norm';
  end
  if ~isfield(opt,'Recenter')
      opt.Recenter = 0;
  end
return
