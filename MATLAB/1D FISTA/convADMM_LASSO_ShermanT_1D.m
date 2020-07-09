function [X_hat, obj, err, l1_norm] = convADMM_LASSO_ShermanT_1D(A0ft_stack,B,X_init,params)
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

[N,M,T] = size(X_init);
BnormSq = sqrt(sum(B.^2, 1));

% B = zeroPad(b,zPad);

% Initialize variables
X_init = forceMaskToZeroArray(X_init,zMask);
Xk = X_init;

% L1 norm variable/lagranage multipliers
Yk = zeros(N,M,T);
Vk = zeros(N,M,T);

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
obj = nan(1,maxIter);
Fit = zeros(N,T);

keep_going = 1;
nIter = 1;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   
    
    % x-update
    for t = 1:T
        Xk(:,:,t) = circulantLinSolve( A0ft_stack,B(:,t),Yk(:,:,t),Vk(:,:,t),params );
    end
    
    % y-update
    Ykp1 = soft(alpha*Xk + (1-alpha)*Yk + Vk,lambda/rho);
    if isNonnegative
        Ykp1(Ykp1<0) = 0;
    end
    % v-update
    Vk = Vk + alpha*Xk + (1-alpha)*Yk - Ykp1;
   
    % Track and display error, objective, sparsity
    for t = 1:T
        Fit(:,t) = Ax_ft_1D(A0ft_stack,Xk);
    end
    
    err(nIter) = 0.5*sum( ((B-Fit).^2)./BnormSq, 'all');
    l1_norm(nIter) = sum(abs(Xk(:)));
    f = err(nIter) + lambda*l1_norm(nIter); %+ rho/2*sum(abs( xkp1(:)-yk(:)+vk(:) ));
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Rho ',     num2str(rho),...
              ' Err ',     num2str(err(nIter)),...
              ' ||x||_1 ', num2str(lambda*l1_norm(nIter)),...
              ' ||x||_0 ', num2str(sum(Xk(:) >0))
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
%         case 'COEF_CHANGE'
%             diff_x = sum(abs(Xkp1(:)-Xk(:)))/numel(Xk);
%             keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end

if adaptRho
    sk = rho*norm(Ykp1(:)-Yk(:));
    rk = norm(Xk(:)-Ykp1(:));
    if rk > mu*sk
        rho = rho*tau;
    elseif sk > mu*rk
        rho = rho/tau;
    end
end

    % Update indices
    Yk = Ykp1;
    if obj(nIter) <= min(obj)
        Xmin = Xk;
    end
end

X_hat = Xmin;
if isNonnegative
    X_hat(X_hat<0) = 0;
end

obj = obj(1:nIter) ;
err = err(1:nIter) ;
l1_norm = l1_norm(1:nIter);

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end
