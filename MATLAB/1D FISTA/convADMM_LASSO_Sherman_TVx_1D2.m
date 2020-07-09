function [x_hat, err, obj, l1_norm, tv_penalty] = convADMM_LASSO_Sherman_TVx_1D2(A0ft_stack,b,x_init,x_n,params)
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

b = zeroPad(b,zPad);
b_ft = fft(b);
bnormsq = sum((b(:)).^2);
tm1 = numel(x_n{1}) > 1;
tp1 = numel(x_n{2}) > 1;

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
xkp1 = x_init;
yk = zeros(size(x_init));
ykp1 = zeros(size(x_init));
z1k = zeros(size(x_init));
z2k = zeros(size(x_init));
vk = zeros(size(xk));
u1k = zeros(size(xk));
u2k = zeros(size(xk));

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
tv_penalty = nan(1,maxIter);
obj = nan(1,maxIter);

% Initial objective
err(1) = sum((b-Ax_ft_1D(A0ft_stack,x_init)).^2);
l1_norm(1) = sum(abs(x_init(:)));
tv_penalty(1) = ( tm1*sum(abs(x_init(:)-x_n{1}(:))) +...
                  tp1*sum(abs(x_init(:)-x_n{2}(:))) );
obj(1) = 0.5/bnormsq*err(1) + lambda*l1_norm(1) + lambda2*tv_penalty(1);

keep_going = 1;
nIter = 1;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   
    
    % x-update
    xkp1 = circulantLinSolveTVx( A0ft_stack,b,b_ft,ykp1,vk,z1k,z2k,u1k,u2k,...
                                 x_n,params);
                
    % y-update
    ykp1 = soft(alpha*xkp1 + (1-alpha)*yk + vk,lambda/(rho));
    if isNonnegative
        ykp1(ykp1<0) = 0;
    end
    % v-update
    vk = vk + alpha*xkp1 + (1-alpha)*yk - ykp1;
    
    % z-update and u-update
    if tm1
        z1kp1 = soft(alpha*xkp1 + (1-alpha)*z1k - x_n{1} + u1k,lambda2/(rho));
        if isNonnegative
            z1kp1(z1kp1<0) = 0;
        end
        u1k = u1k + alpha*xkp1 + (1-alpha)*z1k - z1kp1;
    else
        z1kp1 = 0;
        u1k = 0;
    end
    if tp1
        z2kp1 = soft(alpha*xkp1 + (1-alpha)*z2k - x_n{2} + u2k,lambda2/(rho));
        if isNonnegative
            z2kp1(z2kp1<0) = 0;
        end
        u2k = u2k + alpha*xkp1 + (1-alpha)*z2k - z2kp1;
    else
        z2kp1 = 0;
        u2k = 0;
    end
   
    % Track and display error, objective, sparsity
    fit = Ax_ft_1D(A0ft_stack,xkp1);
        
    err(nIter) = sum((b(:)-fit(:)).^2);
    l1_norm(nIter) = sum(abs(xkp1(:)));
    tv_penalty(nIter) = ( tm1*sum(abs(xkp1(:)-x_n{1}(:))) +...
                          tp1*sum(abs(xkp1(:)-x_n{2}(:))) );
                              
    f = 0.5*err(nIter)/bnormsq + lambda*l1_norm(nIter) + lambda2*tv_penalty(nIter);
        
    obj(nIter) = f;
    if params.verbose
        disp(['Iter ',     num2str(nIter),...
              ' Obj ',     num2str(obj(nIter)),...
              ' Rho ',     num2str(rho),...
              ' Rho2 ',    num2str(rho2),...
              ' Err ',  num2str(0.5*err(nIter)/bnormsq),...
              ' ||x||_1 ', num2str(lambda*l1_norm(nIter)),...
              ' TVx ',     num2str(lambda2*tv_penalty(nIter)),...
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
    if adaptRho
        sk = rho*sum((ykp1(:)-yk(:)).^2);
        skz1 = rho*sum((z1kp1(:)-z1k(:)).^2);
        skz2 = rho*sum((z2kp1(:)-z2k(:)).^2);
        rk = sum((xkp1(:)-ykp1(:)).^2);
        rkz1 = sum((xkp1(:)-z1kp1(:)).^2);
        rkz2 = sum((xkp1(:)-z2kp1(:)).^2);
        if sqrt(rk+rkz1+rkz2) > mu*sqrt(sk+skz1+skz2)
            rho = rho*tau;
        elseif sqrt(sk+skz1+skz2) > mu*sqrt(rk+rkz1+rkz2)
            rho = rho/tau;
        end
    end
%     if adaptRho
%         % rho 1
%         sk = rho*norm(ykp1-yk);
%         rk = norm(xkp1-ykp1);
%         if rk > mu*sk
%             rho = rho*tau;
%         elseif sk > mu*rk
%             rho = rho/tau;
%         end
%         
%         % rho 2
%         skz1 = rho2*norm(z1kp1-z1k);
%         rkz1 = norm(xkp1-z1kp1);
%         skz2 = rho2*norm(z2kp1-z2k);
%         rkz2 = norm(xkp1-z2kp1); 
%         if (tm1&&(rkz1 > mu*skz1)) || (tp1&&(rkz2 > mu*skz2))
%             rho2 = rho2*tau;
%         elseif (tm1&&(skz1 > mu*rkz1)) || (tp1&&(skz2 > mu*rkz2))
%             rho2 = rho2/tau;
%         end
%     end

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
