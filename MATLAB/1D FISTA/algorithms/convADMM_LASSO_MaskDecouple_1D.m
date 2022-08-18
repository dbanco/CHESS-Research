function [x_hat, obj, err, l1_norm] = convADMM_LASSO_MaskDecouple_1D(A0ft_stack,b,x_init,params)
%convADMM_LASSO_1D Image regression by solving LASSO problem 
%                argmin_x 0.5*||Ax-b||^2 + lambda||x||_1
%
% Inputs:
%   b          - (n) polar ring image
%   A0ft_stack - (n x t) fft of unshifted gaussian basis matrices
%   params     - struct containing the following field
%   lambda     - l1 penalty parameter > 0
%   adaptRho   - adaptive rho enable: 1 or 0
%   rho        - admm penalty parameter > 0
%   tau        - adaptive rho parameter: 1.01-1.5
%   mu         - separation factor between primal and dual residual
%   alpha      - momentum parameter1.1-1.8
%   isNonnegative - flag to enforce nonnegative solution
%   x_init - initial guess of solution
%
%   stoppingCriterion - 'OBJECTIVE_VALUE' or 'COEF_CHANGE'
%   tolerance - tolerance for stopping criterion
%   maxIter - maximum number of iterations
%
%   zeroPad         - [row_pad_width,col_pad_width]
%   zeroMask        - [row_indices,col_indices]
%   plotProgress    - 0 or 1
%
% Outputs:
%   x_hat - (n x t) solution 
%   err - (nIters) relative error of solution at each iteration
%   obj - (nIters) objective value of solution at each iteration

tolerance = params.tolerance;
lambda = params.lambda1;
rho = params.rho1;
mu = params.mu;
adaptRho = params.adaptRho;
tau = params.tau;
alpha = params.alpha;
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;

zPad = params.zeroPad;
zMask = params.zeroMask;

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
xkp1 = x_init;

y0k = zeros(size(x_init));
v0k = zeros(size(x_init));

n = numel(b);
y1k = b;
% y1k((1+n-zPad):end) = flipud(b((1+n-2*zPad):(n-zPad)));
% y1k(1:zPad) = flipud(b((1+zPad):(2*zPad)));
y1kp1 = y1k;

v1k = zeros(size(b));

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
obj = nan(1,maxIter);

% Initial objective
err(1) = 0.5*sum((...
            forceMaskToZero(b,zMask)-...
            forceMaskToZero(Ax_ft_1D(A0ft_stack,x_init),zMask) ).^2);
l1_norm(1) = sum(abs(x_init(:)));
obj(1) = err(1) + lambda*l1_norm(1);

keep_going = 1;
nIter = 1;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   
    
    % x-update
    xkp1 = circulantLinSolveMD( A0ft_stack,(y1kp1-v1k),y0k,v0k,params);
%     xkp1 = forceMaskToZeroArray(xkp1,zMask);

    % y-update
    y0kp1 = soft(alpha*xkp1 + (1-alpha)*y0k + v0k,lambda/rho);
    y1kp1 = indepLinSolve(A0ft_stack,b,xkp1,v1k,params);
    if isNonnegative
        y0kp1(y0kp1<0) = 0;
    end
    
    fit = Ax_ft_1D(A0ft_stack,xkp1);
    
    % v-update
    v0k = v0k + alpha*xkp1 + (1-alpha)*y0k - y0kp1;
    v1k = v1k + alpha*fit + (1-alpha)*y1k - y1kp1;
    
    % Track and display error, objective, sparsity
    fitMask = forceMaskToZero(fit,zMask);
    
    err(nIter) = 0.5*sum((forceMaskToZero(b,zMask)-fitMask ).^2);
    l1_norm(nIter) = sum(abs(xkp1(:)));
    f = err(nIter) + lambda*l1_norm(nIter);
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Rho ',     num2str(rho),...
              ' RelErr ',  num2str(err(nIter)),...
              ' ||x||_1 ', num2str(lambda*l1_norm(nIter)),...
              ' ||x||_0 ', num2str(sum(xkp1(:) >0))
               ]);
    end
    
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
    switch params.stoppingCriterion
        case 'OBJECTIVE_VALUE'
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(obj(nIter)-obj(nIter-1));
            keep_going =  (criterionObjective > tolerance);
        case 'COEF_CHANGE'
            diff_x = sum(abs(xkp1(:)-xk(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end

%     if f > prev_f
%         keep_going = 0;
%         xkp1 = xk;
%     end

if adaptRho
    sk = rho*norm(y0kp1-y0k) + rho*norm(y1kp1-y1k);
    rk = norm(xkp1-y0kp1) + norm(y1kp1-fit);
    if rk > mu*sk
        rho = rho*tau;
    elseif sk > mu*rk
        rho = rho/tau;
    end
end

% if criterionObjective < 1e-5
% 	rho = rho/tau;
% end

    % Update indices
    xk = xkp1;
    y0k = y0kp1;
    y1k = y1kp1;
%     v0k = v0kp1;
%     v1k = v1kp1;
    
end

if isNonnegative
    xkp1(xkp1<0) = 0;
end
x_hat = xkp1;
err = err(1:nIter) ;
obj = obj(1:nIter) ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end
