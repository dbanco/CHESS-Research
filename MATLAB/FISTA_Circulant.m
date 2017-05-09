function [x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,b,params)
% Inputs:
% b          - m x n polar ring image
% A0ft_stack - m x n x t x r fft2 of unshifted gaussian basis matrices
% params     - struct containing the following field
%   lambda - l1 penalty parameter > 0
%   L - initial Lipschitz constant > 0
%   beta - backtracking parameter > 1
%   stoppingCriterion - integer indicated stopping criterion (1,2,3)
%   tolerance - tolerance for stopping criterion
%   maxIter - maximum number of iterations
%   isNonnegative - flag to enforce nonnegative solution
%   x_init - initial guess of solution
%
%   Implementation of Fast Iterative Shrinkage-Thresholding Algorithm using 
%   convolutional subroutines for circulant matrix computations as
%   described in:
%   A. Beck and M. Teboulle, “A Fast Iterative Shrinkage-Thresholding 
%       Algorithm for Linear Inverse Problems,” SIAM Journal on Imaging 
%       Sciences, vol. 2, no. 1, pp. 183202, Jan. 2009.

% Define stopping criterion
STOPPING_SPARSE_SUPPORT = 1;
STOPPING_OBJECTIVE_VALUE = 2;
STOPPING_SUBGRADIENT = 3;

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
maxIter = params.maxIter;
isNonnegative = params.isNonnegative;
[m,n,t,r] = size(A0ft_stack) ;
try
    x_init = params.x_init;
catch
    disp('default: x_init = ones(m,n,t,r)')
    x_init = ones(m,n,t,r);
end
if ~all(size(x_init)==[m,n,t,r])
    error('The dimension of the initial xk does not match.');
end

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
l_0 = nan(1,maxIter);

% Initial sparsity and objective
nz_x = (abs(x_init)> eps*10);
f = 0.5*norm(b-Ax_ft_2D(A0ft_stack,x_init))^2 +...
    lambda * norm(x_init(:),1);


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
    grad = AtR_ft_2D(A0ft_stack,Ax_ft_2D(A0ft_stack,ykp1)) - c ; % gradient of f at yk
    
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
        
        % Compute quadratic approximation at yk
        fit2 = Ax_ft_2D(A0ft_stack,ykp1);
        temp1 = norm(b(:)-fit(:))^2;
        temp2 = norm(b(:)-fit2(:))^2 +...
            (xkp1(:)-ykp1(:))'*grad(:) + (L/2)*norm(xkp1(:)-ykp1(:))^2;
        
        % Stop backtrack if objective <= quadratic approximation
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*beta ;
        end
    end
           
    % Track and display error, objective, sparsity
    prev_f = f;
    f = 0.5*norm(b-fit)^2 + lambda * norm(xk(:),1);
    nz_x_prev = nz_x;
    nz_x = (abs(xkp1)>eps*10);
    err(nIter) = norm(b(:)-fit(:))/norm(b(:));
    obj(nIter) = f;
    l_0(nIter) = nz_x;
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
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                keep_going = (criterionActiveSet > tolerance);
            end
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