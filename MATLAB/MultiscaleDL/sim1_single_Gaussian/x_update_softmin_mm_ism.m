function [Xf, cgst, opt] = x_update_softmin_mm_ism(X0, Df, bf, Sf, Y, U, opt, N, M, K, J, T)
% x_update_softmin_reg_gpu -- MM-style linearized softmin solver that uses ISM
%
% Solves the X-subproblem with softmin regularizer by (a) linearizing
% the softmin term at the current X (using compute_softmin to get gradient),
% (b) solving the linear system (A^T A + rho I) X = b_eff exactly with
% your ISM routine `solvedbi_sm`, and (c) repeating for a small number of
% inner linearization (MM) iterations. Optional damping/backtracking
% enforces monotone decrease in the full objective.
%
% USAGE:
%   [Xf, cgst, opt] = x_update_softmin_reg_gpu(X0, Df, bf, Sf, Y, U, opt, N, M, K, J, T)
%
% INPUTS (as in your code)
%   - X0    : initial spatial X (same shape you use, e.g. [1 x N x KJ x T])
%   - Df,bf,Sf,Y,U,opt,... shapes as before
%
% OUTPUT:
%   - Xf    : frequency-domain solution
%   - cgst  : info struct with pit (inner iterations) and last ISM flag if available
%   - opt   : options (unchanged except possibly opt.L or diagnostic fields)
%
% Required functions used:
%   - compute_softmin(X,K,J)  -> [R, dR] where dR is gradient wrt X (spatial)
%   - solvedbi_sm(Df, rho, bf_eff) -> returns Xf that solves (A^T A + rho I) Xf = bf_eff
%       (this is your ISM routine; signature must match)
%

% default options
if ~isfield(opt,'lambda2'), opt.lambda2 = 0; end
if ~isfield(opt,'ism_init'), opt.ism_init = true; end
if ~isfield(opt,'max_inner'), opt.max_inner = 6; end    % number of linearization iterations
if ~isfield(opt,'tol_inner'), opt.tol_inner = 1e-8; end
if ~isfield(opt,'damping'), opt.damping = 1.0; end      % damping factor (1 = full step)
if ~isfield(opt,'backtrack'), opt.backtrack = true; end
if ~isfield(opt,'bt_max'), opt.bt_max = 6; end
if ~isfield(opt,'tau'), opt.tau = 1e-3; end

rho1 = opt.rho1;
lambda2 = opt.lambda2;

% initialize
if isempty(X0)
    % start from zero spatial X if none given
    % dimensions of bf -> use inverse fft to get initial
    try
        X0 = ifft2(bf./(sum(abs(Df).^2,3) + rho1),'symmetric'); % crude
    catch
        X0 = zeros([1, size(bf,1), size(bf,2), size(bf,3)]); % fallback shape - adapt if needed
    end
end

X = X0;                      % spatial domain current iterate
Xf = fft2(X);                % freq domain

% compute initial objective if desired (for backtracking)
compute_obj = @(Xsp) objective_full(Xsp, Df, Sf, Y, U, rho1, lambda2, K, J, opt.tau);

f_cur = compute_obj(X);

% option: initialize with ISM if requested (and lambda2 small)
if opt.ism_init && lambda2 == 0
    % if no softmin, exact ISM directly
    Xf = solvedbi_sm(Df, rho1, bf);
    X = ifft2(Xf, 'symmetric');
    cgst.pit = 1;
    opt.last_inner_flag = 0;
    return;
end

% MM / linearization iterations
cgst.pit = 0;
opt.last_inner_flag = [];

for it = 1:opt.max_inner
    cgst.pit = cgst.pit + 1;
    
    % --- compute gradient of softmin at current X wrt spatial X ---
    % compute_softmin returns [R, dR] where dR is gradient with same shape as X
    [~, dR] = compute_softmin(X, K, J, opt.tau, 'zero', false);
    % ensure dR shape matches X
    if isempty(dR)
        error('compute_softmin did not return gradient dR; required for linearization.');
    end
    
    % build effective RHS in frequency domain:
    % bf is AGSf + rho * fft2(Y - U) (as you use upstream)
    % the linearized normal equations are:
    %   (A^T A + rho I) X = b - lambda2 * dR
    % so bf_eff (freq domain) = bf - lambda2 * fft2(dR)
    bf_eff = bf - lambda2 * fft2(dR);
    
    % --- solve inner linear system exactly with ISM ---
    Xf_new = solvedbi_sm(Df, rho1, bf_eff);
    X_new = ifft2(Xf_new, 'symmetric');
    
    % optional damping/backtracking: accept full step only if objective decreases
    if opt.backtrack
        step = opt.damping;
        f_new = compute_obj(X_new);
        bt = 0;
        while (f_new > f_cur) && (bt < opt.bt_max)
            % damp the update: X_try = X + step*(X_new - X)
            step = step * 0.5;
            X_try = X + step * (X_new - X);
            Xf_try = fft2(X_try);
            % You can reuse ISM to solve for X_try, but here we just evaluate objective on interpolated X.
            f_new = compute_obj(X_try);
            bt = bt + 1;
        end
        % if final step accepted, set X to X_try; else, if even fully damped didn't improve, keep X unchanged
        if f_new <= f_cur
            if bt == 0
                X = X_new; Xf = Xf_new;
            else
                X = X_try; Xf = fft2(X);
            end
            f_cur = f_new;
        else
            % even after damping we couldn't decrease: tighten tolerance or break
            warning('MM+ISM: surrogate solve did not decrease objective after backtracking. Keeping previous X.');
            break;
        end
    else
        % no backtracking: accept full step
        X = X_new; Xf = Xf_new;
        f_cur = compute_obj(X);
    end
    
    opt.last_inner_flag = 0; % if solvedbi_sm returns status, set here
    % stopping criterion (relative reduction)
    if it > 1
        if abs(f_cur - f_prev) < opt.tol_inner * max(1, abs(f_prev))
            break;
        end
    end
    f_prev = f_cur;
end

% final outputs
Xf = fft2(X);
cgst.pit = cgst.pit;
opt.last_obj = f_cur;

end


%% ---------- helper: full objective (spatial) ----------
function f = objective_full(Xsp, Df, Sf, Y, U, rho, lambda2, K, J, tau)
    Xf_local = fft2(Xsp);
    Af = sum(bsxfun(@times, Df, Xf_local), 3);
    data_term = 0.5 * norm(Sf(:) - Af(:))^2;
    rho_term = 0.5 * rho * norm(Xsp(:) - Y(:) + U(:))^2;
    if lambda2 > 0
        R = compute_softmin(Xsp, K, J, tau, 'zero', false);
        if numel(R) > 1
            % compute_softmin may return [R,dR]; handle both outputs
            R = R(1);
        end
    else
        R = 0;
    end
    f = data_term + lambda2 * R + rho_term;
end
