function [x_hat, obj, err, l1_norm] = convADMM_LASSO_Sherman_2D(A0ft_stack,b,x_init,params)
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

b = zeroPad(b,zPad);
bnormsq = sum((b(:)).^2);

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
xkp1 = x_init;
yk = zeros(size(x_init));
ykp1 = zeros(size(x_init));
vk = zeros(size(xk));

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
obj = nan(1,maxIter);

% Initial objective
err(1) = sum((b-Ax_ft_2D(A0ft_stack,x_init)).^2,'all');
l1_norm(1) = sum(abs(x_init(:)));
obj(1) = 0.5/bnormsq*err(1) + lambda*l1_norm(1);

keep_going = 1;
nIter = 1;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   
    
    % x-update
    xkp1 = circulantLinSolve_2D( A0ft_stack,b,ykp1,vk,params );
%     xkp1 = forceMaskToZeroArray(xkp1,zMask);

    % y-update
    ykp1 = soft(alpha*xkp1 + (1-alpha)*yk + vk,lambda/rho);
    if isNonnegative
        ykp1(ykp1<0) = 0;
    end
    % v-update
    vk = vk + alpha*xkp1 + (1-alpha)*yk - ykp1;
   
    % Track and display error, objective, sparsity
    fit = Ax_ft_2D(A0ft_stack,xkp1);
    
    err(nIter) = sum((b(:)-fit(:)).^2)/(2*bnormsq);
    l1_norm(nIter) = sum(abs(xkp1(:)));
    f = err(nIter) + lambda*l1_norm(nIter);
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Rho ',     num2str(rho),...
              ' Err ',  num2str(err(nIter)),...
              ' ||x||_1 ', num2str(lambda*l1_norm(nIter)),...
              ' ||x||_0 ', num2str(sum(xkp1(:) >0))
               ]);
    end
    
    if params.plotProgress
        figure(1)    
        subplot(2,1,1)
        imagesc(b)
        subplot(2,1,2)
        imagesc(fit)
        
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
    sk = rho*norm(ykp1(:)-yk(:));
    rk = norm(xkp1(:)-ykp1(:));
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
    yk = ykp1;
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
