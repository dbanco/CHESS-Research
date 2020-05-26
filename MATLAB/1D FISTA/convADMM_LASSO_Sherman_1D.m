function [x_hat, err, obj] = convADMM_LASSO_Sherman_1D(A0ft_stack,b,x_init,params)
%convADMM_LASSO_1D Image regression by solving LASSO problem 
%                argmin_x 0.5*||Ax-b||^2 + lambda||x||_1
%
% Inputs:
% b          - (n) polar ring image
% A0ft_stack - (n x t) fft of unshifted gaussian basis matrices
% params     - struct containing the following field
%   lambda - l1 penalty parameter > 0
%   rho - admm penalty parameter > 0
%   L - gradient step size > 0
%   isNonnegative - flag to enforce nonnegative solution
%   x_init - initial guess of solution

%   stoppingCriterion - integer indicated stopping criterion (1,2,3)
%   tolerance - tolerance for stopping criterion
%   maxIter - maximum number of iterations
%   maxGradIter - maximum number of iterations for gradient descent

%
% Outputs:
% x_hat - (n x t) solution 
% err - (nIters) relative error of solution at each iteration
% obj - (nIters) objective value of solution at each iteration
% l_0 - (nIters) sparsity of solution at each iteration


% Define stopping criterion
STOPPING_OBJECTIVE_VALUE = 1;
COEF_CHANGE = 2;

% Get parameters
stoppingCriterion = params.stoppingCriterion;
tolerance = params.tolerance;
lambda = params.lambda;
rho = params.rho;
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;

zPad = params.zeroPad;
zMask = params.zeroMask;
[n,t] = size(A0ft_stack) ;

if ~all(size(x_init)==[n,t])
    error('The dimension of the x_init does not match.');
end

b = zeroPad(b,zPad);
bnorm = norm(b(:));

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
obj = nan(1,maxIter);


% Initial objective
f = 0.5/bnorm*sum((b-Ax_ft_1D(A0ft_stack,x_init)).^2) +...
    lambda * sum(abs(x_init(:)));
prev_f = f;

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
xkp1 = x_init;
yk = x_init;
vk = zeros(size(xk));

keep_going = 1;
nIter = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   
    
    % x-update
    xkp1 = circulantLinSolve( A0ft_stack,b,yk,params );

    % nonegativity
    if isNonnegative
    xkp1(xkp1<0) = 0;
    end

    % y-update
    yk = soft(xkp1 + vk,lambda/rho);
    
    % v-update
    vk = vk + xkp1 - yk;
   
    % Track and display error, objective, sparsity
    fit = Ax_ft_1D(A0ft_stack,xkp1);
    
    err(nIter) = sum((b(:)-fit(:)).^2)/2/bnorm;
    l1_norm(nIter) = lambda*sum(abs(xkp1(:)));
    f = err(nIter) + l1_norm(nIter); %+ rho/2*sum(abs( xkp1(:)-yk(:)+vk(:) ));
        
    obj(nIter) = f;
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' RelErr ',  num2str(err(nIter)),...
          ' ||x||_1 ', num2str(l1_norm(nIter)),...
          ' ||x||_0 ', num2str(sum(xkp1(:) >0))
           ]);
    
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
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f-prev_f);
            keep_going =  (criterionObjective > tolerance);
        case COEF_CHANGE
            diff_x = sum(abs(xkp1(:)-xk(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end

%     if f > prev_f
%         keep_going = 0;
%         xkp1 = xk;
%     end

    % Update indices
    prev_f = f;
    xk = xkp1; 
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
