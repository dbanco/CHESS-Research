function [x_hat, err, obj] = convADMM_LASSO_Sherman_TVx_1D(A0ft_stack,b,x_init,x_n,params)
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
mu = params.mu;
adaptRho = params.adaptRho;
tau = params.tau;
alpha = params.alpha;
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;

rho2 = params.rho2;
lambda2 = params.lambda2;

zPad = params.zeroPad;
zMask = params.zeroMask;
[n,t] = size(A0ft_stack) ;

if ~all(size(x_init)==[n,t])
    error('The dimension of the x_init does not match.');
end
if numel(x_n) == 1
    xEnd = 1;
else
    xEnd = 0;
end
b = zeroPad(b,zPad);
bnormsq = sum((b(:)).^2);

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
xkp1 = x_init;
yk = x_init;
ykp1 = x_init;
z1k = x_init;
z1kp1 = x_init;
z2k = zeros(size(x_init));
z2kp1 = x_init;
vk = zeros(size(xk));
u1k = zeros(size(xk));
u2k = zeros(size(xk));

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
tv_penalty = nan(1,maxIter);
obj = nan(1,maxIter);

% Initial objective
err(1) = 0.5/bnormsq*sum((b-Ax_ft_1D(A0ft_stack,x_init)).^2);
l1_norm(1) = lambda*sum(abs(x_init(:)));
obj(1) = err(1) + l1_norm(1);

keep_going = 1;
nIter = 1;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   
    
    % x-update
    xkp1 = circulantLinSolveTVx( A0ft_stack,b,ykp1,vk,z1k,z2k,u1k,u2k,params,xEnd );

    % y-update
    ykp1 = soft(alpha*xkp1 + (1-alpha)*yk + vk,lambda/rho);
    if isNonnegative
        ykp1(ykp1<0) = 0;
    end
    
    % z-update
    if xEnd
        x_n{2} = xkp1;
    end
    z1kp1 = soft(xkp1 - x_n{1} + u1k,lambda2/rho2);
    z2kp1 = soft(xkp1 - x_n{2} + u2k,lambda2/rho2);
    if isNonnegative
        z1kp1(z1kp1<0) = 0;
    end
    if isNonnegative
        z2kp1(z2kp1<0) = 0;
    end
    
    % v-update
    vk = vk + alpha*xkp1 + (1-alpha)*yk - ykp1;
    
    % u-update
    u1k = u1k + xkp1 - z1kp1;
    u2k = u2k + xkp1 - z2kp1;
   
    % Track and display error, objective, sparsity
    fit = Ax_ft_1D(A0ft_stack,xkp1);
        
    err(nIter) = sum((b(:)-fit(:)).^2)/2/bnormsq;
    l1_norm(nIter) = lambda*sum(abs(xkp1(:)));
    tv_penalty(nIter) = lambda2*( sum(abs(xkp1(:)-x_n{1}(:))) +...
                                  sum(abs(xkp1(:)-x_n{2}(:))) );
                              
    f = err(nIter) + l1_norm(nIter) + tv_penalty(nIter);
        
    obj(nIter) = f;
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' Rho ',     num2str(rho),...
          ' Rho2 ',    num2str(rho2),...
          ' RelErr ',  num2str(err(nIter)),...
          ' ||x||_1 ', num2str(l1_norm(nIter)),...
          ' TVx ',     num2str(tv_penalty(nIter)),...
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
            criterionObjective = abs(obj(nIter)-obj(nIter-1));
            keep_going =  (criterionObjective > tolerance);
        case COEF_CHANGE
            diff_x = sum(abs(xkp1(:)-xk(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end

    if adaptRho
        % rho 1
        sk = rho*norm(ykp1-yk);
        rk = norm(xkp1-ykp1);
        if rk > mu*sk
            rho = rho*tau;
        elseif sk > mu*rk
            rho = rho/tau;
        end
        
        % rho 2
        skz1 = rho2*norm(z1kp1-z1k);
        rkz1 = norm(xkp1-z1kp1);
        if rkz1 > mu*skz1
            rho2 = rho2*tau;
        elseif skz1 > mu*rkz1
            rho2 = rho2/tau;
        end
        skz2 = rho2*norm(z2kp1-z2k);
        rkz2 = norm(xkp1-z2kp1);
        if rkz2 > mu*skz2
            rho2 = rho2*tau;
        elseif skz2 > mu*rkz2
            rho2 = rho2/tau;
        end   
    end

    % Update indices
    xk = xkp1;
    yk = ykp1;
    z1k = z1kp1;
    z2k = z2kp1;
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
