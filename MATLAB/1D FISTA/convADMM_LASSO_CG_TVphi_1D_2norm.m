function [X_hat, err, obj, l1_norm, l2_penalty] = convADMM_LASSO_CG_TVphi_1D_2norm(A0ft_stack,B,X_init,params)
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
%

% Get parameters
tolerance = params.tolerance;
lambda1 = params.lambda1;
lambda2 = params.lambda2;
rho1 = params.rho1;
rho2 = params.rho2;
mu = params.mu;
adaptRho = params.adaptRho;
tau = params.tau;
alpha = params.alpha;
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;

zPad = params.zeroPad;
zMask = params.zeroMask;

[N,K,T] = size(X_init);
BnormSq = sum(B.^2, 1);

% Initialize variables
X_init = forceMaskToZeroArray(X_init,zMask);
Xk = X_init;

% L1 norm variable/lagranage multipliers
Yk = zeros(N,K,T);
Vk = zeros(N,K,T);

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
l2_penalty = nan(1,maxIter);
obj = nan(1,maxIter);

keep_going = 1;
nIter = 0;
count = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1;   
    
    % x-update
    [Xkp1,cgIters] = conjGrad_TVphi_1D_2norm( A0ft_stack,B,Xk,(Yk-Vk),params);

    % y-update and v-update
    for t = 1:T
        Ykp1(:,:,t) = soft(alpha*params.lambda1(t)*Xkp1(:,:,t) + (1-alpha)*Yk(:,:,t) + Vk(:,:,t), lambda1(t)/rho1);
    end
    if isNonnegative
        Ykp1(Ykp1<0) = 0;
    end
    for t = 1:T
        Vk(:,:,t) = Vk(:,:,t) + alpha*params.lambda1(t)*Xkp1(:,:,t) + (1-alpha)*Yk(:,:,t) - Ykp1(:,:,t);
    end
    % Track and display error, objective, sparsity
    fit = Ax_ft_1D_Time(A0ft_stack,Xkp1);
    err(nIter) = sum(((B-fit).^2)./BnormSq,'all');
    Xsum = 0;
    for t = 1:T
        Xsum = Xsum + lambda1(t)*sum(abs(Xkp1(:,:,t)),'all');
    end
    l1_norm(nIter) = Xsum;
    l2_penalty(nIter) = lambda2*sum((DiffPhiX_1D(Xkp1).^2),'all');
    f = 0.5*err(nIter) + l1_norm(nIter) + l2_penalty(nIter);
    
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' cgIters ',  num2str(cgIters),...
              ' Rho1 ',     num2str(rho1),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Err ',     num2str(0.5*err(nIter)),...
              ' ||x||_1 ', num2str(l1_norm(nIter)),...
              ' TVx ',     num2str(l2_penalty(nIter)),...
              ' ||x||_0 ', num2str(sum(Xkp1(:) >0))
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
            diff_x = sum(abs(Xkp1(:)-Xk(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end
    
    if adaptRho
        skY = rho1*sum((Ykp1(:)-Yk(:)).^2);
        rkY = sum((Xkp1(:)-Ykp1(:)).^2);
        if sqrt(rkY) > mu*sqrt(skY)
            rho1 = rho1*tau;
        elseif sqrt(skY) > mu*sqrt(rkY)
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

%     if obj(nIter) > 1e10
%         keep_going = 0;
%     end

    % Update indices
    Xk = Xkp1;
    Yk = Ykp1;
    if obj(nIter) <= min(obj)
        Xmin = Xkp1;
    end
    
end

X_hat = Xmin;
if isNonnegative
    X_hat(X_hat<0) = 0;
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
