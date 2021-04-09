function [x_hat, err, obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_1D(A0ft_stack,b,x_init,params)
%convADMM_LASSO_1D Image regression by solving LASSO problem 
%                argmin_x 0.5*||Ax-b||^2 + lambda||x||_1
%
% Inputs:
%   b          - (n) polar ring image
%   A0ft_stack - (n x t) fft of unshifted gaussian basis matrices
%   params     - struct containing the following field
%   lambda     - l1 penalty parameter > 0
%   lambda2    - TVx penalty parameter > 0
%   adaptRho   - adaptive rho enable: 1 or 0
%   rho        - admm penalty parameter > 0
%   rho2       - admm penalty parameter > 0
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
% x_hat - (n x t) solution 
% err - (nIters) relative error of solution at each iteration
% obj - (nIters) objective value of solution at each iteration
% l_0 - (nIters) sparsity of solution at each iteration

% Get parameters
tolerance = params.tolerance;
lambda1 = params.lambda1;
rho1 = params.rho1;
mu = params.mu;
adaptRho = params.adaptRho;
tau = params.tau;
alpha = params.alpha;
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;

zPad = params.zeroPad;
zMask = params.zeroMask;

[N,M] = size(x_init);

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;

% L1 norm variable/lagranage multipliers
yk = zeros(N,M);
ykp1 = zeros(N,M);
vk = zeros(N,M);

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
tv_penalty = nan(1,maxIter);
obj = nan(1,maxIter);

keep_going = 1;
nIter = 0;
count = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1;   
    
    % x-update
    [xkp1,cgIters] = conjGrad_1D( A0ft_stack,b,xk,(yk-vk),params);
    
    % y-update and v-update
    ykp1 = soft(alpha*xkp1 + (1-alpha)*yk + vk, lambda1/(rho1));
    if isNonnegative
        ykp1(ykp1<0) = 0;
    end
    vk = vk + alpha*xkp1 + (1-alpha)*yk - ykp1;
    
    % Track and display error, objective, sparsity
    fit = Ax_ft_1D(A0ft_stack,xkp1);
    
    err(nIter) = sum(((b-fit).^2),'all');
    l1_norm(nIter) = sum(abs(xkp1),'all');
    
    f = 0.5*err(nIter) + lambda1*l1_norm(nIter);
    
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' cgIters ', num2str(cgIters),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Rho1 ',     num2str(rho1),...
              ' Err ',     num2str(0.5*err(nIter)),...
              ' ||x||_1 ', num2str(lambda1*l1_norm(nIter)),...
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
            try
                criterionObjective = abs(obj(nIter)-obj(nIter-1));
                keep_going =  (criterionObjective > tolerance);
            catch
                keep_going = 1;
            end
        case 'COEF_CHANGE'
            diff_x = sum(abs(xkp1(:)-xk(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end

    if adaptRho
        skY = rho1*norm(ykp1(:)-yk(:));
        rkY =norm(xkp1(:)-ykp1(:));
        if rkY > mu*skY
            rho1 = rho1*tau;
        elseif skY > mu*rkY
            rho1 = rho1/tau;
        end
    end
    
    % Stop if objective increases consistently
    if (nIter > 10) && (obj(nIter-1) < obj(nIter))
        count = count + 1;
        if count > 20
            keep_going = 0;
        end
    else
        count = 0;
    end
    
    % Update indices
    xk = xkp1;
    yk = ykp1;
    if obj(nIter) <= min(obj)
        xMin = xkp1;
    end
    
end

x_hat = xMin;
if isNonnegative
    x_hat(x_hat<0) = 0;
end

err = err(1:nIter) ;
obj = obj(1:nIter) ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end
