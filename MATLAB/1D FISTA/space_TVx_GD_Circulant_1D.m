function [x_hat, err, t_k, L, obj, l_0] = space_TVx_GD_Circulant_1D(A0ft_stack,b,neighbors_x,x_init,params)
%FISTA_Circulant Image regression by solving LASSO problem 
%                argmin_x ||Ax-b||^2 + lambda||x|| +...
%                         gamma sum_{adjacent_xi}^2 (1/4)||xn-x||^2
%
%   Implementation of Fast Iterative Shrinkage-Thresholding Algorithm using 
%   convolutional subroutines for circulant matrix computations as
%   described in:
%   A. Beck and M. Teboulle, “A Fast Iterative Shrinkage-Thresholding 
%       Algorithm for Linear Inverse Problems,�? SIAM Journal on Imaging 
%       Sciences, vol. 2, no. 1, pp. 183202, Jan. 2009.
%
% Inputs:
% b          - (n) polar ring images
% neighibors_ev - (scalar) average expected variance of neighboring ring images
% A0ft_stack - (n x t) fft2 of unshifted gaussian basis matrices
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
% x_hat - (n x t) solution 
% err - (nIters) relative error of solution at each iteration
% obj - (nIters) objective value of solution at each iteration
% l_0 - (nIters) sparsity of solution at each iteration


% Define stopping criterion
STOPPING_OBJECTIVE_VALUE = 1;
COEF_CHANGE = 2;
GRADIENT_NORM = 3;

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
tvBeta = params.tvBeta;
beta = params.beta;
maxIter = params.maxIterReg;
isNonnegative = params.isNonnegative;

zPad = params.zeroPad;
zMask = params.zeroMask;

tvBeta = params.tvBeta;
numIms = params.numIms;
imageNum = params.imageNum;

[n,t] = size(A0ft_stack) ;
if ~all(size(x_init)==[n,t])
    error('The dimension of the initial xk does not match.');
end

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
obj1 = nan(1,maxIter);
obj2 = nan(1,maxIter);
obj3 = nan(1,maxIter);
l_0 = nan(1,maxIter);

b = zeroPad(b,zPad);
bnorm = norm(b);

if numel(neighbors_x) == 2
    D =  [-1 1 0;... 
          0 -1 1];
elseif imageNum == 1
    D = [-1 1];
elseif imageNum == numIms
    D = [1 -1];
end

% Initial sparsity and objective
f_obj = 0.5/bnorm*norm(b-forceMaskToZero(Ax_ft_1D(A0ft_stack,x_init),zMask))^2 +...
        lambda*sum(abs(x_init(:)));

% Add reg term  
tvObj = 0;
for i = 1:numel(neighbors_x)
    tvObj = tvObj + sum(sqrt( (x_init(:)-neighbors_x{i}(:)).^2 + tvBeta^2 ));
end
f_obj = f_obj + params.gamma*tvObj;

prev_L = 0.5;
prev_f = f_obj;
min_f = f_obj;
min_x = x_init;
min_iter = 0;
old_count = 0;

% Used to compute gradient
c = AtR_ft_1D(A0ft_stack,b)/bnorm;

x_init = forceMaskToZeroArray(x_init,zMask);
xkm1 = x_init;
xk = x_init;
zk = xk;
t_k = params.t_k;
t_kp1 = params.t_k;
keep_going = 1;
nIter = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;
    
    % Data matching gradient update
    grad = AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,zk),zMask))/bnorm...
        - c + lambda*zk./sqrt(zk.^2 + tvBeta^2);;

    % TV regularizer gradient update    
    if numel(neighbors_x) == 2
    	f_tv = [neighbors_x{1}(:), zk(:),  neighbors_x{2}(:)]'; 
    elseif imageNum == 1
        f_tv = [zk(:), neighbors_x{1}(:)]'; 
    elseif imageNum == numIms
        f_tv = [zk(:), neighbors_x{1}(:)]';
    end
    
    gradTV = params.gamma*gradientTV(f_tv,tvBeta,D);
    grad = grad + reshape(gradTV,size(zk));
    

    % Backtracking
    stop_backtrack = 0 ;
    while ~stop_backtrack 
        xk = zk - (1/L)*grad ;
        if isNonnegative
            xk(xk<0) = 0;
        end
        
        % Compute objective at xk
        fit_xk = forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask);
        tvObj_xk = 0;
        for i = 1:numel(neighbors_x)
            tvObj_xk = tvObj_xk + sum(sqrt( (xk(:)-neighbors_x{i}(:)).^2 + tvBeta^2 ));
        end
        
        f_xk = 0.5/bnorm*norm(b(:)-fit_xk(:))^2 +...
              lambda*sum(sqrt(xk(:).^2 + tvBeta^2)) +...
              params.gamma*tvObj_xk;
        
        % Compute objective at zk
        fit_zk = forceMaskToZero(Ax_ft_1D(A0ft_stack,zk),zMask); 
        tvObj_zk = 0;
        for i = 1:numel(neighbors_x)
            tvObj_zk = tvObj_zk + sum(sqrt( (zk(:)-neighbors_x{i}(:)).^2 + tvBeta^2 ));
        end
        
        f_zk = 0.5/bnorm*norm(b(:)-fit_zk(:))^2 +...
            lambda*sum(sqrt(zk(:).^2 + tvBeta^2)) +...
            params.gamma*tvObj_zk;
        
        criterion = f_zk - (0.5/L)*sum(grad(:).^2); 
        
        % Stop backtrack if objective <= quadratic approximation
        if params.noBacktrack
            stop_backtrack = 1 ;
        elseif f_xk <= criterion
            stop_backtrack = 1;
            params.noBacktrack = 1;
        else
            L = L*beta ;
            if L > 1e15
                keep_going = 0;
                stop_backtrack = 1;
            end
        end
    end
    
    % Track and display error, objective, sparsity
    prev_f = f_obj;
    f_data = 0.5/bnorm*norm(b-fit_xk)^2;
    f_sparse = lambda*sum(abs(xk(:)));
%     f_sparse = lambda*sum(sqrt(xk(:).^2 + tvBeta^2));
    f_tv = params.gamma*tvObj_xk;
    f_obj = f_data + f_sparse + f_tv;   
    err(nIter) = norm(b(:)-fit_xk(:));
    obj(nIter) = f_obj;
    obj1(nIter) = f_data;
    obj2(nIter) = f_sparse;
    obj3(nIter) = f_tv;
    l_0(nIter) = sum(abs(xk(:))>0);
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' L ',       num2str(L),...
          ' ||x||_0 ', num2str(l_0(nIter)),...
          ' ||x||_1 ',  num2str(sum(abs(xk(:)))),...
          ' RelErr ',  num2str(err(nIter)),...
          ' Obj1 ',     num2str(obj1(nIter)),...
          ' Obj2 ',     num2str(obj2(nIter)),...
          ' Obj3 ',     num2str(obj3(nIter))]);
    
	if params.plotProgress
        lim1 = 0;
        lim2 = max(b(:));
        figure(1)
       
        subplot(2,3,1)
        plot(b);
        title('img')
        
        subplot(2,3,2)
        plot(Ax_ft_1D(A0ft_stack,xk));
        title('xk')
        
        subplot(2,3,3)
        plot(fit_zk);
        title('zk')
        
        subplot(2,3,4)
        fit_gk = forceMaskToZero(Ax_ft_1D(A0ft_stack,gk),zMask);
        plot(fit_gk);
        title('gk')
        
        subplot(2,3,5)
        plot(abs(b-fit));
        title('diff xk')
        
        subplot(2,3,6)
        plot(abs(b-fit_zk));
        title('diff zk')
        
        pause(0.05)
    end
      
    % Check stopping criterion
    switch stoppingCriterion
        case GRADIENT_NORM
            keep_going = norm(grad(:)) > tolerance;
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f_obj-prev_f);
            keep_going =  (criterionObjective > tolerance);
        case COEF_CHANGE
            diff_x = sum(abs(xk(:)-xkm1(:)))/numel(xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end
    if f_obj < min_f
        min_x = xk;
        min_f = f_obj;
        min_iter = nIter;
        old_count = 0;
    else
        old_count = old_count + 1;
        if old_count > 50
            keep_going = 0;
        end
    end
    % Update indices
    prev_f = f_obj;
    zk = xk;
end

x_hat = min_x;
err = err(1:min_iter) ;
obj = obj(1:min_iter) ;
l_0 = l_0(1:min_iter) ;