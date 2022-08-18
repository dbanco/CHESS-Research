function [x_hat err obj l_0] = GD_test(A,b,x_init,params)
%FISTA_Circulant_1D Image regression by solving LASSO problem 
%                argmin_x 0.5*||Ax-b||^2 + lambda||x||_1
%
%   Implementation of Fast Iterative Shrinkage-Thresholding Algorithm using 
%   convolutional subroutines for circulant matrix computations as
%   described in:
%   A. Beck and M. Teboulle, “A Fast Iterative Shrinkage-Thresholding 
%       Algorithm for Linear Inverse Problems,�? SIAM Journal on Imaging 
%       Sciences, vol. 2, no. 1, pp. 183202, Jan. 2009.
%
% Inputs:
% b          - (n) polar ring image
% A0ft_stack - (n x t) fft of unshifted gaussian basis matrices
% params     - struct containing the following field
%   lambda - l1 penalty parameter > 0
%   L - initial Lipschitz constant > 0
%   beta - backtracking parameter > 1
%   stoppingCriterion - integer indicated stopping criterion (1,2,3)
%   tolerance - tolerance for stopping criterion
%   maxIter - maximum number of iterations
%   isNonnegative - flag to enforce nonnegative solution
%   x_init - initial guess of solution
%
% Outputs:
% x_hat - (n x t) solution 
% err - (nIters) relative error of solution at each iteration
% obj - (nIters) objective value of solution at each iteration
% l_0 - (nIters) sparsity of solution at each iteration


% Define stopping criterion
GRADIENT_NORM = 1;
STOPPING_OBJECTIVE_VALUE = 2;
COEF_CHANGE = 3;

% Set default parameter values
% stoppingCriterion = STOPPING_OBJECTIVE_VALUE;
% maxIter = 200 ;
% isNonnegative = 1;
% tolerance = 1e-3;
% x_init = ones(m,n,t,r) ;

% Get parameters
stoppingCriterion = params.stoppingCriterion;
tolerance = params.tolerance;
L = params.L;
lambda = params.lambda;
beta = params.beta;
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;

tvBeta = params.tvBeta;
numIms = params.numIms;
imageNum = params.imageNum;

zPad = params.zeroPad;
zMask = params.zeroMask;
[n,t] = size(A) ;

b = zeroPad(b,zPad);
bnorm = norm(b(:));

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
l_0 = nan(1,maxIter);

% Initial sparsity and objective
f = 0.5*sum((b-A*x_init).^2)+...
    lambda * sum(sqrt(x_init(:).^2 + tvBeta^2));

prev_L = 0.5;
prev_f = f;
min_f = f;
old_count = 0;


x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
zk = xk;
keep_going = 1;
nIter = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;        
    
    % Compute gradient of f
    grad = A'*(A*xk - b) + lambda*xk./sqrt(xk.^2 + tvBeta^2);
    
    stop_backtrack = 0;
    while ~stop_backtrack 
        
        %l1/nonnegative-proximal
        zk = xk - (1/L)*grad ;
%         imagesc(zk)
        if isNonnegative
            zk(zk<0) = 0;
        end
        
        % Compute objective at xk/zk and criterion to stop backtracking
        f_xk = 0.5*sum((b-A*xk).^2) + lambda*sum(sqrt(xk(:).^2 + tvBeta^2));
        fit = A*zk;
        f = 0.5*sum((b-fit).^2) + lambda*sum(sqrt(zk(:).^2 + tvBeta^2));
        
        criterion = f_xk - (0.25/L)*sum(grad(:).^2);

        
        % Stop backtrack logic
        if params.noBacktrack
            stop_backtrack = 1 ;
        elseif f <= criterion
            stop_backtrack = 1;
            params.noBacktrack = 1;
        else
            L = L*beta ;
            if L > 1e50
                keep_going = 0;
                stop_backtrack = 1;
            end
        end
    end
   
    % Track and display error, objective, sparsity
    err(nIter) = norm(b(:)-fit(:));
    obj(nIter) = f;
    l_0(nIter) = sum(abs(zk(:))>eps*10);
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' L ',       num2str(L),...
          ' ||x||_0 ', num2str(l_0(nIter)),...
          ' ||x||_1 ',  num2str(sum(abs(zk(:)))),...
          ' RelErr ',  num2str(err(nIter)) ]);
    
    if params.plotProgress

        figure(1)    
        hold off
        plot(b)
        hold on
        plot(fit)
        legend('data','fit')
        
        pause
    end
    
    % Check stopping criterion
    switch stoppingCriterion
        case GRADIENT_NORM
            keep_going = norm(grad(:)) > tolerance;
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f-prev_f);
            keep_going =  (criterionObjective > tolerance);
        case COEF_CHANGE
            diff_x = sum(abs(zk(:)-xk(:)))/numel(zk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end
    if f < min_f
        min_x = zk;
        min_f = f;
        min_iter = nIter;
        old_count = 0;
    else
        old_count = old_count + 1;
        if old_count > 20
            keep_going = 0;
        end
    end
    
    % Update indices
    prev_f = f;
    xk = zk;
    
end

x_hat = min_x;
err = err(1:min_iter) ;
obj = obj(1:min_iter) ;
l_0 = l_0(1:min_iter) ;