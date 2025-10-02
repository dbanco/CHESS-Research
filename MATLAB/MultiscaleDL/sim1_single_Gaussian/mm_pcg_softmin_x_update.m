function [Xf, cgst, opt] = mm_pcg_softmin_x_update(X0, Df, bf, Sf, Y, U, opt, N, M, K, J, T)
% MM_PCG_SOFTMIN_X_UPDATE
%   MM (strict: linearize phi in r) + PCG solve of surrogate quadratic
%
% Usage:
%   [Xf, cgst, opt] = mm_pcg_softmin_x_update(X0, Df, bf, Sf, Y, U, opt, N, M, K, J, T)
%
% Inputs:
%   X0    : initial spatial X (shape [1, N, KJ, T]) (can be empty)
%   Df    : frequency-domain dictionary [Nx Ny M] with M = K*J
%   bf    : RHS in frequency domain (AGSf + rho*fft2(Y-U)) [Nx Ny M]
%   Sf    : S in frequency domain (or spatial S depending on your objective eval)
%   Y,U   : spatial arrays same shape as X
%   opt   : options struct with fields used below
%   N,M,K,J,T : dims (M should equal K*J)
%
% Outputs:
%   Xf    : frequency domain solution
%   cgst  : struct with fields {flag, iter, relres, mm_iters}
%   opt   : options updated with last objective etc.
%
% Required subfunctions inside: compute_softmin_weights, reg_op_spatial, objective_eval
%
% Important options (opt):
%   opt.rho1      - rho
%   opt.lambda2   - softmin weight
%   opt.tau       - softmin temperature
%   opt.mm_maxit  - number of MM iterations (default 6)
%   opt.mm_tol    - tolerance on surrogate objective change
%   opt.pcgTol    - PCG tolerance (default 1e-8)
%   opt.pcgMaxIter- PCG max iterations
%   opt.precond_eps - small eps added to preconditioner diagonal
%   opt.backtrack - whether to require objective decrease and re-solve with tighter tol
%
% Notes:
%   - This code uses a frequency-domain diagonal preconditioner approximating
%     P(omega,m) = |Df(omega,m)|^2 + rho + lambda2 * P_reg(omega,m).
%   - The computation of P_reg uses the shift masks and the MM weights and maps
%     to channels m = (k-1)*J + j. If your channel mapping differs, adjust accordingly.

% ----------------------
% Default options
if ~isfield(opt,'rho1'), error('opt.rho1 required'); end
if ~isfield(opt,'lambda2'), opt.lambda2 = 0; end
if ~isfield(opt,'tau'), opt.tau = 1e-3; end
if ~isfield(opt,'mm_maxit'), opt.mm_maxit = 6; end
if ~isfield(opt,'mm_tol'), opt.mm_tol = 1e-8; end
if ~isfield(opt,'pcgTol'), opt.pcgTol = 1e-8; end
if ~isfield(opt,'pcgMaxIter'), opt.pcgMaxIter = 2000; end
if ~isfield(opt,'precond_eps'), opt.precond_eps = 1e-8; end
if ~isfield(opt,'backtrack'), opt.backtrack = true; end
% ----------------------

rho1 = opt.rho1;
lambda2 = opt.lambda2;
tau = opt.tau;

[Nx, Ny, Mdf] = size(Df);
if Mdf ~= M
    warning('Df third dim (%d) != M (%d). Proceeding with M=%d from Df.', Mdf, M);
    M = Mdf;
end
vec = @(x) x(:);

% init X (spatial) and Xf
if isempty(X0)
    % crude initialization: if bf available, simple inverse diag
    P0 = sum(abs(Df).^2, 3) + rho1; % [Nx Ny]
    try
        Xf0 = bsxfun(@rdivide, bf, P0);
        X0 = ifft2(Xf0,'symmetric');
    catch
        X0 = zeros(1, N, K*J, T);
    end
end
X = X0;
Xf = fft2(X);

% prepare Aop and Atop (freq domain)
Aop = @(Xf_in) sum(bsxfun(@times, Df, Xf_in), 3);  % returns [Nx Ny] complex
Atop = @(rf) ifft2(bsxfun(@times, conj(Df), rf),'symmetric'); % rf [Nx Ny] -> returns [Nx Ny M] spatial after ifft

% objective evaluator
objective_eval_local = @(Xsp) objective_eval(Xsp, Df, Sf, Y, U, rho1, lambda2, K, J, tau);

% compute bf_spatial once
b_spatial = ifft2(bf,'symmetric');  % spatial domain of bf

% MM iterations
cgst = struct('flag',[],'iter',[],'relres',[],'mm_iters',0);
f_old = objective_eval_local(X);

% create R_op action using helper reg_op_spatial
function y = L_op(u_vec)
    Usp = unvec_sp(u_vec);
    % data term
    Ufx = fft2(Usp);
    AUf = Aop(Ufx);  % [Nx Ny]
    AtAU = ifft2(bsxfun(@times, conj(Df), AUf),'symmetric'); % spatial per-channel: same as Atop(AUf)
    % AtAU computed as sum(conj(Df).*AUf) per channel then ifft2; implement via Atop(AUf)
    % add rho term
    Yrho = rho1 * Usp;
    % add reg term (linear operator for fixed weights)
    Rterm = lambda2 * reg_op_spatial(Usp, W_all, K, J, 'zero', false);
    Ysp = AtAU + Yrho + Rterm;
    y = vec_sp(Ysp);
end

% preconditioner function M^{-1} acting on vector v (spatial domain)
function z = Mfun(v)
    Vsp = unvec_sp(v);
    Vf = fft2(Vsp);  % frequency
    Zf = bsxfun(@rdivide, Vf, P_approx); % per-channel division
    Zsp = ifft2(Zf,'symmetric');
    z = vec_sp(Zsp);
end

for mm = 1:opt.mm_maxit
    cgst.mm_iters = cgst.mm_iters + 1;
    % 1) compute softmin weights W at current X
    [W_all, Rsum] = compute_softmin_weights(X, K, J, tau, 'zero', false);
    % W_all: K x L x Nsp x (T-1) in our helper; but we will average over pixels for precond.
    % We'll use reg_op_spatial with W_all for the operator action.
    
    % 2) build frequency-domain approximate diagonal for surrogate reg term: P_reg(omega,m)
    %    we compute per-channel P_reg by summing |1 - exp(-i omega·shift)|^2 times
    %    aggregated weights for the corresponding atom k.
    P_reg = compute_reg_diag_freq(W_all, K, J, Nx, Ny); % returns [Nx Ny M]
    
    % 3) build preconditioner diagonal approx: P_approx = |Df|^2 + rho + lambda2 * P_reg
    P_approx = abs(Df).^2 + rho1 + lambda2 * P_reg;
    P_approx = P_approx + opt.precond_eps; % numerical safety
    Pvec = P_approx(:);  % flattened for Mfun
    
    % 4) build PCG operator in spatial domain: L(u) = Atop(Aop(fft2(u))) + rho*u + lambda2 * Rop(u; W_all)
    xsz = size(X);
    nvec = numel(X);
    unvec_sp = @(v) reshape(v, xsz);
    vec_sp = @(Z) reshape(Z, [], 1);

    % RHS is b_spatial (spatial domain)
    rhs_vec = vec_sp(b_spatial);
    x0_vec = vec_sp(X); % warm start
    
    % run PCG
    tol = opt.pcgTol;
    maxit = opt.pcgMaxIter;
    [xv, flag, relres, iter] = pcg(@(u)L_op(u), rhs_vec, tol, maxit, @(v)Mfun(v), [], x0_vec);
    
    X_sol = reshape(xv, xsz);
    Xf_new = fft2(X_sol);
    
    % Evaluate objective, require decrease (optional) - re-run with tighter PCG tol if not decreased
    f_new = objective_eval_local(X_sol);
    if opt.backtrack && (f_new > f_old + 1e-12)
        % re-run with stricter tolerance
        tol2 = max(tol*1e-3, 1e-12);
        maxit2 = 2*maxit;
        [xv2, flag2, relres2, iter2] = pcg(@(u)L_op(u), rhs_vec, tol2, maxit2, @(v)Mfun(v), [], x0_vec);
        X_sol2 = reshape(xv2, xsz);
        f_new2 = objective_eval_local(X_sol2);
        if f_new2 <= f_old + 1e-12
            X_sol = X_sol2; Xf_new = fft2(X_sol); f_new = f_new2;
            flag = flag2; relres = relres2; iter = iter2;
        else
            % if still worse, keep previous X (reject) and warn
            warning('MM-PCG surrogate solve did not decrease objective; rejecting update this MM iteration.');
            X_sol = X; Xf_new = fft2(X_sol); f_new = f_old;
        end
    end
    
    % accept
    X = X_sol;
    Xf = Xf_new;
    f_old = f_new;
    cgst.flag = flag; cgst.iter = iter; cgst.relres = relres;
    
    % termination check
    if mm > 1
        if abs(f_old - f_prev) < opt.mm_tol * max(1, abs(f_prev))
            break;
        end
    end
    f_prev = f_old;
end

% outputs
opt.last_obj = f_old;
cgst.mm_iters = mm;
Xf = fft2(X);

end