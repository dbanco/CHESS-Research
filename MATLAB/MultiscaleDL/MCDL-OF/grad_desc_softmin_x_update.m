function [Xf, objVals, L] = grad_desc_softmin_x_update(Xf0,Df,Sf,Y,U,opt,K,J,M,opts)
% Gradient descent solver for
%   1/2||b - A x||^2 + softmin(x) + (rho1/2)||x - z + w||^2
%
% Inputs:
%   b, A      : data and operator
%   z, w      : ADMM variables
%   rho1      : quadratic penalty weight
%   x0        : initial x
%   grad_softmin : function handle @(x) -> gradient of softmin term
%   opts      : struct with fields:
%               maxIter, tol, L0 (initial step), backtracking
%
% Output:
%   x         : solution
%   objVals   : objective history

if ~isfield(opts,'maxIter'), opts.maxIter = 300; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'L0'), opts.L0 = 1e1; end
if ~isfield(opts,'backtracking'), opts.backtracking = true; end
if ~isfield(opts,'Lmin'), opts.Lmin = 1e-12; end
if ~isfield(opts,'Lmax'), opts.Lmax = 1e12; end


Xf = Xf0;
X = ifft2(Xf,'symmetric');

objVals = zeros(opts.maxIter,1);
rho1 = opt.rho1;
lambda2 = opt.lambda2;

L_spec = max(abs(Df(:)).^2) + rho1;
L = max(L_spec, opt.L);    % start conservatively at spectral

Aop = @(xf) sum(bsxfun(@times,Df,xf),3);
Atop = @(rf) ifft2(bsxfun(@times, conj(Df), rf),'symmetric');

obj_softmin = @(x) compute_softmin(x, K, J);

f_cur = objective(Xf,X,Sf,Aop,Y,U,rho1,lambda2,obj_softmin);
objVals(1) = f_cur;

for k = 1:opts.maxIter
    % --- Gradient of smooth part ---
    
    rf = Aop(Xf) - Sf;
    % [~,gradSoftMin] = compute_softmin(X,K,J);
    % Clip gradients of softmin
    % gradSoftMin = min(max(gradSoftMin, -1e0), 1e0);
    grad = Atop(rf) + rho1*(X - Y + U); %+ lambda2*gradSoftMin
    
    step = 1/L;
    
    % --- Backtracking ---
    if opts.backtracking
        f_prev = f_cur;
        max_backtracks = 60;
        bt = 0;
        while true
            X_new = X - step*grad;
            Xf_new = fft2(X_new);
            f_new = objective(Xf_new,X_new,Sf,Aop,Y,U,rho1,lambda2,obj_softmin);
            % Quadratic upper bound
            diff = X_new(:)-X(:);
            f_quad = f_prev + grad(:).'*diff + (L/2)*(diff'*diff);
            if ((f_new <= f_quad) || bt >= max_backtracks) 
                break;
            else
                L = min(2*L, opts.Lmax);
                step = 1/L;
                bt = bt + 1;
            end
        end
    else
        X_new = X - step*grad;
        f_new = objective(Xf_new,X_new,Sf,Aop,Y,U,rho1,lambda2,obj_softmin);
    end

    % fprintf('iter %d: f=%.6e, ||grad||=%.3e, L=%.3e\n',k,f_cur,norm(grad(:)),L);
    
    if f_new <= f_cur
        % successful step: small reduction (speed up)
        L = max(L/2, opts.Lmin);
    end

    % --- Update ---
    denom = max(norm(X(:)), 1e-12);
    rel_change = norm(X_new(:)-X(:))/denom;
 
    if rel_change < opts.tol
        objVals = objVals(1:k);
        break;
    end
    
    if bt < max_backtracks
        X = X_new;
        Xf = fft2(X);
    end
    f_cur = f_new;
    objVals(k) = f_cur;
end

end

%% Objective function
function f = objective(Xf,X,Sf,Aop,Y,U,rho1,lambda2,softmin_func)
    f = 0.5*norm(Sf - Aop(Xf),'fro')^2;% + lambda2*softmin_func(X); 
    f = f + (rho1/2)*norm(X - Y + U,'fro')^2;
end

