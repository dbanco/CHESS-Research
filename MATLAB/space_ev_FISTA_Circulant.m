function [x_hat, err, obj, l_0] = space_ev_FISTA_Circulant(A0ft_stack,b,neighbors_ev,var_theta,x_init,params)
%FISTA_Circulant Image regression by solving LASSO problem 
%                argmin_x ||Ax-b||^2 + lambda||x|| +...
%                         gamma sum_{adjacent_xi}^4 (1/4)||xn-x||^2
%
%   Implementation of Fast Iterative Shrinkage-Thresholding Algorithm using 
%   convolutional subroutines for circulant matrix computations as
%   described in:
%   A. Beck and M. Teboulle, “A Fast Iterative Shrinkage-Thresholding 
%       Algorithm for Linear Inverse Problems,�? SIAM Journal on Imaging 
%       Sciences, vol. 2, no. 1, pp. 183202, Jan. 2009.
%
% Inputs:
% b          - (m x n) polar ring images
% neighibors_ev - (scalar) average expected variance of neighboring ring images
% A0ft_stack - (m x n x t x r) fft2 of unshifted gaussian basis matrices
% params     - struct containing the following field
%   lambda - l1 penalty parameter > 0
%   gamma - spatial smoothness regularization parameter > 0
%   L - initial Lipschitz constant > 0
%   beta - backtracking parameter > 1
%   stoppingCriterion - integer indicated stopping criterion (1,2,3)
%   tolerance - tolerance for stopping criterion
%   maxIter - maximum number of iterations
%   isNonnegative - flag to enforce nonnegative solution
%   x_init - initial guess of solution
%
% Outputs:
% x_hat - (m x n x t x r) solution 
% err - (nIters) relative error of solution at each iteration
% obj - (nIters) objective value of solution at each iteration
% l_0 - (nIters) sparsity of solution at each iteration


% Define stopping criterion
STOPPING_OBJECTIVE_VALUE = 1;
STOPPING_SUBGRADIENT = 2;

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
params.gamma = params.gamma*sum(b(:));
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;
[m,n,t,r] = size(A0ft_stack) ;

if ~all(size(x_init)==[m,n,t,r])
    error('The dimension of the initial xk does not match.');
end

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
l_0 = nan(1,maxIter);

% Initial sparsity and objective
f = 0.5*norm(b-Ax_ft_2D(A0ft_stack,x_init))^2 +...
    lambda * norm(x_init(:),1);

% Add spatial expected variance smoothness part of objective
x_ev = compute_exp_az_variance(x_init,var_theta);
f = f + params.gamma*(x_ev - neighbors_ev)^2;

% Used to compute gradient
c = AtR_ft_2D(A0ft_stack,b);

xkm1 = x_init;
xk = x_init;
t_k = 1;
t_kp1 = 1;
keep_going = 1;
nIter = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;
    
    ykp1 = xk + ((t_k-1)/t_kp1)*(xk-xkm1) ;
    if isNonnegative
        ykp1(ykp1<0) = 0;
    end
        
    % Compute gradient of f
    x_ev = compute_exp_az_variance(ykp1,var_theta);
    total = sum(ykp1(:));
    grad = AtR_ft_2D(A0ft_stack,Ax_ft_2D(A0ft_stack,ykp1)) - c;
    for az_i = 1:numel(var_theta)
    	grad(:,:,az_i,:) = grad(:,:,az_i,:) + params.gamma*(x_ev - neighbors_ev)*((var_theta(az_i)-x_ev)/total); 
    end
    % Backtracking
    stop_backtrack = 0 ;
    while ~stop_backtrack 
        gk = ykp1 - (1/L)*grad ;
        if isNonnegative
            gk(gk<0) = 0;
        end
        xkp1 = soft(gk,lambda/L) ;
        
        % Compute objective at xkp1
        fit = Ax_ft_2D(A0ft_stack,xkp1);
        x_ev = compute_exp_az_variance(xkp1,var_theta);
        temp1 = norm(b(:)-fit(:))^2 + params.gamma*(x_ev - neighbors_ev)^2;
        
        % Compute quadratic approximation at ykp1
        fit2 = Ax_ft_2D(A0ft_stack,ykp1);
        x_ev = compute_exp_az_variance(ykp1,var_theta);
        temp2 = norm(b(:)-fit2(:))^2 +...
            params.gamma*(x_ev - neighbors_ev)^2 +...
            (xkp1(:)-ykp1(:))'*grad(:) +...
            (L/2)*norm(xkp1(:)-ykp1(:))^2;
        
        % Stop backtrack if objective <= quadratic approximation
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*beta ;
        end
    end
           
    % Track and display error, objective, sparsity
    prev_f = f;
    f = temp1 + lambda * norm(xk(:),1);
    err(nIter) = norm(b(:)-fit(:))/norm(b(:));
    obj(nIter) = f;
    l_0(nIter) = sum(abs(xkp1(:))>eps*10);
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' ||x||_0 ', num2str(l_0(nIter)),...
          ' RelErr ',  num2str(err(nIter)) ]);
    
    % Check stopping criterion
    switch stoppingCriterion
        case STOPPING_SUBGRADIENT
            sk = L*(ykp1-xkp1) +...
                 AtR_ft_2D(A0ft_stack,Ax_ft_2D(A0ft_stack,(xkp1-ykp1)));
            keep_going = norm(sk(:)) > tolerance*L*max(1,norm(xkp1(:)));
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f-prev_f)/(prev_f);
            keep_going =  (criterionObjective > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end
    
    t_k = t_kp1;
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k));
    xkm1 = xk;
    xk = xkp1;
    
    
end

x_hat = xkp1 ;
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