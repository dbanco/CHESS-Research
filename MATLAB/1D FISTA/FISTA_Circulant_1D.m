function [x_hat, err, obj, l_0, t_k, L] = FISTA_Circulant_1D(A0ft_stack,b,x_init,params)
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
STOPPING_OBJECTIVE_VALUE = 1;
STOPPING_SUBGRADIENT = 2;
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
zPad = params.zeroPad;
zMask = params.zeroMask;
[n,t] = size(A0ft_stack) ;

if ~all(size(x_init)==[n,t])
    error('The dimension of the initial x_init does not match.');
end

b = zeroPad(b,zPad);
bnormsq = sum(b(:).^2);

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
l_0 = nan(1,maxIter);

% Initial sparsity and objective
f = 0.5/bnormsq*norm(b-Ax_ft_1D(A0ft_stack,x_init))^2 +...
    lambda * norm(x_init(:),1);

% Used to compute gradient
c = AtR_ft_1D(A0ft_stack,b)/bnormsq;

x_init = forceMaskToZeroArray(x_init,zMask);
xkm1 = x_init;
xk = x_init;
zk = xk;
t_k = 1;
t_kp1 = 1;
keep_going = 1;
nIter = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;        
    
    % Compute gradient of f
    grad = AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,zk),zMask))/bnormsq - c; % gradient of f at zk
    
    stop_backtrack = 0;
    while ~stop_backtrack 
        
        %l1/nonnegative-proximal
        gk = zk - (1/L)*grad ;
        xk = soft(gk,lambda/L) ;
        if isNonnegative
            xk(xk<0) = 0;
        end
        
        % Compute objective at xk
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask);
        
        % Compute quadratic approximation at yk
        fit2 = forceMaskToZero(Ax_ft_1D(A0ft_stack,zk),zMask);
        temp1 = 0.5*sum((b(:)-fit(:)).^2)/bnormsq  + lambda*sum(abs(xk(:)));
        temp2 = 0.5*sum((b(:)-fit2(:)).^2)/bnormsq + lambda*sum(abs(zk(:))) +...
            (xk(:)-zk(:))'*grad(:) + (L/2)*norm(xk(:)-zk(:))^2;
        
        % Stop backtrack if objective <= quadratic approximation
        if params.noBacktrack
            stop_backtrack = 1 ;
        elseif temp1 <= temp2
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
    
    
    t_kp1 = 0.5*(1 + sqrt(1+4*t_k^2));
    zk = xk + ((t_k-1)/t_kp1)*(xk-xkm1);    

    % Track and display error, objective, sparsity
    prev_f = f;
    f = 0.5/bnormsq*norm(b-fit)^2 + lambda * sum(abs(xk(:)));
    err(nIter) = norm(b(:)-fit(:))^2;
    obj(nIter) = f;
    l_0(nIter) = sum(abs(zk(:))>eps*10);
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' L ',       num2str(L),...
          ' ||x||_0 ', num2str(l_0(nIter)),...
          ' ||x||_1 ',  num2str(sum(abs(xk(:)))),...
          ' Err ',  num2str(err(nIter)) ]);
    
    if params.plotProgress
        lim1 = 0;
        lim2 = max(b(:));
        figure(1)
       
        subplot(2,3,1)
        plot(b);
        title('img')
        
        subplot(2,3,2)
        plot(Ax_ft_1D(A0ft_stack,xk));
        title('xk')
        
        subplot(2,3,3)
        plot(fit2);
        title('zk')
        
        subplot(2,3,4)
        fit_gk = forceMaskToZero(Ax_ft_1D(A0ft_stack,gk),zMask);
        plot(fit_gk);
        title('gk')
        
        subplot(2,3,5)
        plot(abs(b-fit));
        title('diff xk')
        
        subplot(2,3,6)
        plot(abs(b-fit2));
        title('diff zk')
        
        pause(0.05)
    end
    

    % Check stopping criterion
    switch stoppingCriterion
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f-prev_f);
            keep_going =  (criterionObjective > tolerance);
        case STOPPING_SUBGRADIENT
            sk = L*(xk-xkm1) +...
                 AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,xk-xkm1),zMask))/bnorm;
            keep_going = norm(sk(:)) > tolerance*L*max(1,norm(xk(:)));
        case COEF_CHANGE
            diff_x = sum(abs(xk(:)-xkm1(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end
    
    % Update indices
    t_k = t_kp1;
    xkm1 = xk;
end

x_hat = xk;
err = err(1:nIter) ;
obj = obj(1:nIter) ;
l_0 = l_0(1:nIter) ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end