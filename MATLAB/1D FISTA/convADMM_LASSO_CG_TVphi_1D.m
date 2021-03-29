function [X_hat, err, obj, l1_norm, tv_penalty] = convADMM_LASSO_CG_TVphi_1D(A0ft_stack,B,X_init,params)
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
Bnorms = zeros(T,1);
for j = 1:T
    Bnorms(j) = norm(B(:,j));
    B(:,j) = B(:,j)/Bnorms(j);
end

% Initialize variables
X_init = forceMaskToZeroArray(X_init,zMask);
Xk = X_init;

% L1 norm variable/lagranage multipliers
Yk = zeros(N,K,T);
Ykp1 = zeros(N,K,T);
Vk = zeros(N,K,T);

% TVx variable/lagranage multipliers
Zk = zeros(K,T-1);
Uk = zeros(K,T-1);

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
    [Xkp1,cgIters] = conjGrad_TVphi_1D( A0ft_stack,B,Bnorms,Xk,(Yk-Vk),(Zk-Uk),params,zMask);
    
    % y-update and v-update
    for t = 1:T
        Ykp1(:,:,t) = soft(alpha*Xkp1(:,:,t) + (1-alpha)*Yk(:,:,t) + Vk(:,:,t), lambda1(t)/rho1);
    end
    if isNonnegative
        Ykp1(Ykp1<0) = 0;
    end
    Vk = Vk + alpha*Xkp1 + (1-alpha)*Yk - Ykp1;
    
    % z-update and u-update
    Zkp1 = soft(DiffPhiX_1D(Xkp1) + Uk, lambda2/rho2);
    Uk = Uk + DiffPhiX_1D(Xkp1) - Zkp1;
    
    % Track and display error, objective, sparsity
    fit = Ax_ft_1D_Time(A0ft_stack,Xkp1,Bnorms);
    err(nIter) = sum(((B-fit).^2),'all');
    Xsum = 0;
    for t = 1:T
        Xsum = Xsum + lambda1(t)*sum(abs(Xkp1(:,:,t)),'all');
    end
    l1_norm(nIter) = Xsum;
    tv_penalty(nIter) = lambda2*sum(abs(DiffPhiX_1D(Xkp1)),'all');
    f = 0.5*err(nIter) + l1_norm(nIter) + tv_penalty(nIter);
    
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' cgIters ',  num2str(cgIters),...
              ' Rho1 ',     num2str(rho1),...
              ' Rho2 ',     num2str(rho2),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Err ',     num2str(0.5*err(nIter)),...
              ' ||x||_1 ', num2str(l1_norm(nIter)),...
              ' TVx ',     num2str(tv_penalty(nIter)),...
              ' ||x||_0 ', num2str(sum(Xkp1(:) >0))
               ]);
    end
    
    if params.plotProgress
        figure(1)    
        hold off
        plot(B(:,10))
        hold on
        plot(fit(:,10))
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
        skZ = rho2*sum((Zkp1(:)-Zk(:)).^2);
        rkY = sum((Xkp1(:)-Ykp1(:)).^2);
        rkZ = sum((DiffPhiX_1D(Xkp1) - Zkp1).^2,'all');
        if sqrt(rkY) > mu*sqrt(skY)
            rho1 = rho1*tau;
        elseif sqrt(skY) > mu*sqrt(rkY)
            rho1 = rho1/tau;
        end
        if sqrt(rkZ) > mu*sqrt(skZ)
            rho2 = rho2*tau;
        elseif sqrt(skZ) > mu*sqrt(rkZ)
            rho2 = rho2/tau;
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
    Zk = Zkp1;
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
