function [D, Y, Gh, Gv, optinf] = cbpdndl_sep(D0, D0h, D0v, S, lambda, opt)

%-------------------------------------------------------------------------------------------------------%
%   Author: Jorge Quesada <jorge.quesada@pucp.edu.pe>                         			        %
%   This code is based on the SPORCO library framework (http://brendt.wohlberg.net/software/SPORCO/),   %
%   and thus makes use of several subroutines in said library, as well as maintains its structure and   %
%   variable nomenclature.										%
%-------------------------------------------------------------------------------------------------------%

% cbpdndl_sep -- Convolutional BPDN Separable Dictionary Learning
%
%         argmin_{x_m,dh_m, dv_m} (1/2) \sum_k ||\sum_m dh_m * dv_m * x_k,m - s_k||_2^2 +
%                           lambda \sum_k \sum_m ||x_k,m||_1
%
%         The solution is computed using Augmented Lagrangian methods
%         (see boyd-2010-distributed) with efficient solution of the 
%         main linear systems (see wohlberg-2016-efficient).
%
% Usage:
%       [D, Y, Gh, Gv, optinf] = cbpdndl_sep(D0, D0h, D0v, S, lambda, opt)
%
% Input:
%       D0          Initial dictionary
%       D0v         Vertical component of initial dictionary
%       D0h         Horizontal component of initial dictionary
%       S           Input images
%       lambda      Regularization parameter
%       opt         Options/algorithm parameters structure (see below)
%
% Output:
%       D           Dictionary filter set (3D array)
%       Y           Coefficient maps (4D array)
%	Gh          Horizontal components of dictionary filter set
%       Gv          Vertical components of dictionary filter set
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


if nargin < 4,
  opt = [];
end
checkopt(opt, defaultopts([]));
opt = defaultopts(opt);

% Set up status display for verbose operation
hstr = ['Itn   Fnc       DFid      l1        Cnstr     '...
        'r(X)      s(X)      r(D)      s(D) '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e';
nsep = 84;
if opt.AutoRho,
  hstr = [hstr '     rho  '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
if opt.AutoSigma,
  hstr = [hstr '     sigma  '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
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
if size(S,3) > 1,
  xsz = [size(S,1) size(S,2) size(D0,3) size(S,3)];
  % Insert singleton 3rd dimension (for number of filters) so that
  % 4th dimension is number of images in input s volume
  S = reshape(S, [size(S,1) size(S,2) 1 size(S,3)]);
else
  xsz = [size(S,1) size(S,2) size(D0,3) 1];
end
Nx = prod(xsz);
Nd = prod(xsz(1:2))*size(D0,3);
cgt = opt.CGTol;

% Dictionary size may be specified when learning multiscale
% dictionary
if isempty(opt.DictFilterSizes),
  dsz = [size(D0,1) size(D0,2)];
else
  dsz = opt.DictFilterSizes;
end

% Mean removal and normalisation projections
Pzmn = @(x) bsxfun(@minus, x, mean(mean(x,1),2));
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));
Pnrmh = @(x) bsxfun(@rdivide, x, sqrt(mean(sum(x.^2, 2), 1)));
Pnrmv = @(x) bsxfun(@rdivide, x, sqrt(mean(sum(x.^2, 1), 2)));

% Projection of filter to full image size and its transpose
% (zero-pad and crop respectively)
Pzp = @(x) zpad(x, xsz(1:2));
PzpT = @(x) bcrop(x, dsz);

Pzpv = @(x) zpad(x, [xsz(1) size(x,2)]);
Pzph = @(x) zpad(x, [size(x,1) xsz(2)]);

PzpvT = @(x) bcrop(x, [dsz(1) size(x,2)]);
PzphT = @(x) bcrop(x, [size(x,1) dsz(2)]);
% Projection of dictionary filters onto constraint set
if opt.ZeroMean,
  Pcn = @(x) Pnrm(Pzp(Pzmn(PzpT(x))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
  Pcnv = @(x) Pnrmv(Pzpv(PzpvT(x)));
  Pcnh = @(x) Pnrmh(Pzph(PzphT(x)));
end

% Start timer
tstart = tic;

% Project initial dictionary onto constraint set
D = Pnrm(D0);

% Compute signal in DFT domain
Sf = fft2(S);
Sfv = fft(S);
Sfh = fft(S,[],2);
% Set up algorithm parameters and initialise variables
rho = opt.rho;
if isempty(rho), rho = 50*lambda+1; end;
if opt.AutoRho,
  asgr = opt.RhoRsdlRatio;
  asgm = opt.RhoScaling;
end
sigmah = opt.sigma; sigmav = opt.sigma; sigma=opt.sigma;
if isempty(sigmah), sigmah = size(S,3); end;
if opt.AutoSigma,
  asdr = opt.SigmaRsdlRatio;
  asdm = opt.SigmaScaling;
end
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
X = [];
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
Df = [];
if isempty(opt.G0),
  G = Pzp(D);
else
  G = opt.G0;
end
Gprv = G;
if isempty(opt.H0),
  if isempty(opt.G0),
    H = zeros(size(G), class(S));
  else
    H = G;
  end
else
  H = opt.H0;
end

Gf = fft2(G, size(S,1), size(S,2));
GSf = bsxfun(@times, conj(Gf), Sf);


Gh = Pzph(D0h);
Ghprv=Gh;

Gfh = fft(Gh, [], 2);
Dfh = Gfh;
Hh = zeros(size(Gh), class(S));

Gv = Pzpv(D0v);
Gvprv=Gv;

Hv = zeros(size(Gv), class(S));
Gfv = fft(Gv);
Dfv = Gfv;

% Main loop
k = 1;
while k <= opt.MaxMainIter && (rx > eprix|sx > eduax|rd > eprid|sd >eduad),

  % Solve X subproblem. It would be simpler and more efficient (since the
  % DFT is already available) to solve for X using the main dictionary
  % variable D as the dictionary, but this appears to be unstable. Instead,
  % use the projected dictionary variable G
  Xf = solvedbi_sm(Gf, rho, GSf + rho*fft2(Y - U));
  X = ifft2(Xf, 'symmetric');
  clear Xf Gf GSf;

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
  Yfh = fft(Y, [],2);
  Yfv = fft(Y);
  
  %YSf = sum(bsxfun(@times, conj(Yf), Sf), 4);

  % Update dual variable corresponding to X, Y
  U = U + Xr - Y;
  clear Xr;

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
  clear X;

  % Compute l1 norm of Y
  Jl1 = sum(abs(vec(bsxfun(@times, opt.L1Weight, Y))));

  % Update record of previous step Y
  Yprv = Y;


  % Solve Dh and Dv subproblems (1 admm pass for each side). Similarly, it would be simpler and more efficient to
  % solve for D using the main coefficient variable X as the coefficients,
  % but it appears to be more stable to use the shrunk coefficient variable Y
  %----Dv----%
  Yaux = ifft(bsxfun(@times, Gfh, Yfh), [], 2, 'symmetric');
  YSf = sum(sum(bsxfun(@times, conj(fft(Yaux)), Sfv), 2),4);
  
  [Dfv, cgst] = solvemdbi_cgv(fft(Yaux), sigmav, YSf + sigmav*fft(Gv - Hv), ...
                              cgt, opt.MaxCGIter, Dfv(:));
 
  clear YSf;
  Dv = ifft(Dfv, 'symmetric');

  % See pg. 21 of boyd-2010-distributed
  if opt.DRelaxParam == 1,
    Drv = Dv;
  else
    Drv = opt.DRelaxParam*Dv + (1-opt.DRelaxParam)*Gv;
  end

  % Solve Gv subproblem
  Gv = Pcnv(Drv + Hv);
  Gfv = fft(Gv);

  % Update dual variable corresponding to Dv, Gv
  Hv = Hv + Drv - Gv;
  clear Drv;

  %----Dh----%
  Yaux = ifft(bsxfun(@times, Gfv, Yfv), 'symmetric');
  YSf = sum(sum(bsxfun(@times, conj(fft(Yaux, [], 2)), Sfh), 1),4);

  [Dfh, cgst] = solvemdbi_cgh(fft(Yaux,[],2), sigmah, YSf + sigmah*fft(Gh - Hh,[],2), ...
                              cgt, opt.MaxCGIter, Dfh(:));

  clear YSf;
  Dh = ifft(Dfh,[],2, 'symmetric');



  % See pg. 21 of boyd-2010-distributed
  if opt.DRelaxParam == 1,
    Drh = Dh;
  else
    Drh = opt.DRelaxParam*Dh + (1-opt.DRelaxParam)*Gh;
  end

  % Solve Gh subproblem
  Gh = Pcnh(Drh + Hh);
  Gfh = fft(Gh, [], 2);


  % Update dual variable corresponding to Dh, Gh
  Hh = Hh + Drh - Gh;
  clear Drh;

  for c=1:size(Gh,3)
      G(:,:,c) = Gv(:,1,c)*Gh(1,:,c);
      Dr(:,:,c) = Dv(:,:,c)*Dh(:,:,c);
  end
  Gf = fft2(Pcn(G));
  GSf = bsxfun(@times, conj(Gf), Sf);
  D=Dr; 
  
  clear Dr;
  % Compute primal and dual residuals and stopping thresholds for D update
  nD = norm(D(:)); nG = norm(G(:)); nH = norm(H(:));
  nDh = norm(Dh(:)); nGh = norm(Gh(:)); nHh = norm(Hh(:));
  nDv = norm(Dv(:)); nGv = norm(Gv(:)); nHv = norm(Hv(:));
  if opt.StdResiduals,
    % See pp. 19-20 of boyd-2010-distributed
    rd = norm(vec(D - G));
    sd = norm(vec(sigma*(Gprv - G)));
    eprid = sqrt(Nd)*opt.AbsStopTol+max(nD,nG)*opt.RelStopTol;
    eduad = sqrt(Nd)*opt.AbsStopTol+sigma*nH*opt.RelStopTol;
  else
    % See wohlberg-2015-adaptive
    rd = norm(vec(D - G))/max(nD,nG);
    sd = norm(vec(Gprv - G))/nH;
    eprid = sqrt(Nd)*opt.AbsStopTol/max(nD,nG)+opt.RelStopTol;
    eduad = sqrt(Nd)*opt.AbsStopTol/(sigma*nH)+opt.RelStopTol;
    

    rdh = norm(vec(Dh - Gh))/max(nDh,nGh);
    sdh = norm(vec(Ghprv - Gh))/nHh;
    %eprid = sqrt(Nd)*opt.AbsStopTol/max(nDh,nGh)+opt.RelStopTol;
    %eduad = sqrt(Nd)*opt.AbsStopTol/(sigmah*nHh)+opt.RelStopTol;

    
    rdv = norm(vec(Dv - Gv))/max(nDv,nGv);
    sdv = norm(vec(Gvprv - Gv))/nHv;
    %eprid = sqrt(Nd)*opt.AbsStopTol/max(nDh,nGh)+opt.RelStopTol;
    %eduad = sqrt(Nd)*opt.AbsStopTol/(sigmah*nHh)+opt.RelStopTol;
  end

  % Apply CG auto tolerance policy if enabled
  if opt.CGTolAuto && (rd/opt.CGTolFactor) < cgt,
    cgt = rd/opt.CGTolFactor;
  end

  % Compute measure of D constraint violation
  Jcn = norm(vec(Pcn(D) - D));
  clear D;

  % Update record of previous step G
  Gprv = G;
  Gvprv=Gv; Ghprv= Gh;

  % Compute data fidelity term in Fourier domain (note normalisation)
  Jdf = sum(vec(abs(sum(bsxfun(@times,Gf,Yf),3)-Sf).^2))/(2*xsz(1)*xsz(2));
  clear Yf;
  Jfn = Jdf + lambda*Jl1;


  % Record and display iteration details
  tk = toc(tstart);
  optinf.itstat = [optinf.itstat;...
       [k Jfn Jdf Jl1 rx sx rd sd eprix eduax eprid eduad rho sigma tk]];
  if opt.Verbose,
    dvc = [k, Jfn, Jdf, Jl1, Jcn, rx, sx, rd, sd];
    if opt.AutoRho,
      dvc = [dvc rho];
    end
    if opt.AutoSigma,
      dvc = [dvc sigma];
    end
    disp(sprintf(sfms, dvc));
  end

  % See wohlberg-2015-adaptive and pp. 20-21 of boyd-2010-distributed
  if opt.AutoRho,
    if k ~= 1 && mod(k, opt.AutoRhoPeriod) == 0,
      if opt.AutoRhoScaling,
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
  if opt.AutoSigma,
    if k ~= 1 && mod(k, opt.AutoSigmaPeriod) == 0,
      if opt.AutoSigmaScaling,
        sigmlth = sqrt(rdh/sdh);
        sigmltv = sqrt(rdv/sdv);
        if sigmlth < 1, sigmlth = 1/sigmlth; end
        if sigmlth > opt.SigmaScaling, sigmlth = opt.SigmaScaling; end
               
        if sigmltv < 1, sigmltv = 1/sigmltv; end
        if sigmltv > opt.SigmaScaling, sigmltv = opt.SigmaScaling; end
      else
        sigmlth = opt.SigmaScaling; sigmltv = opt.SigmaScaling;
      end
      ssfh = 1; ssfv = 1;
      if rdh > opt.SigmaRsdlRatio*sdh, ssfh = sigmlth; end
      if sdh > opt.SigmaRsdlRatio*rdh, ssfh = 1/sigmlth; end
      
      if rdv > opt.SigmaRsdlRatio*sdv, ssfv = sigmltv; end
      if sdv > opt.SigmaRsdlRatio*rdv, ssfv = 1/sigmltv; end
      sigmah = ssfh*sigmah;  sigmav = ssfv*sigmav;
      Hh = Hh/ssfh; Hv = Hv/ssfv;
    end
  end


  k = k + 1;

end

D = PzpT(G);

% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Y = Y;
optinf.U = U;
optinf.G = G;
optinf.H = H;
optinf.lambda = lambda;
optinf.rho = rho;
optinf.sigma = sigma;
optinf.cgt = cgt;
if exist('cgst'), optinf.cgst = cgst; end

if opt.Verbose && opt.MaxMainIter > 0,
  disp(char('-' * ones(1,nsep)));
end
Gh=PzphT(Gh(:,:,:));
Gv=PzpvT(Gv(:,:,:));
Y=Y;
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


function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'AbsStopTol'),
    opt.AbsStopTol = 1e-6;
  end
  if ~isfield(opt,'RelStopTol'),
    opt.RelStopTol = 1e-4;
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

return
