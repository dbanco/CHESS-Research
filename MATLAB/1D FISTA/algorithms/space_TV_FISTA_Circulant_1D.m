function [x_hat, err, t_k, L, obj, l_0] = space_TV_FISTA_Circulant_1D(A0ft_stack,b,neighbors_vdf,x_init,params)
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
STOPPING_SUBGRADIENT = 2;
COEF_CHANGE = 3;

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

if numel(neighbors_vdf) == 2
    D =  [-1 1 0;... 
          0 -1 1];
elseif imageNum == 1
    D = [-1 1];
elseif imageNum == numIms
    D = [1 -1];
end

% Initial sparsity and objective
f_obj = 0.5/bnorm*norm(b-forceMaskToZero(Ax_ft_1D(A0ft_stack,x_init),zMask))^2 +...
    lambda * norm(x_init(:),1);

% Add vdf tv reg term  
vdf = squeeze(sum(x_init,1));
vdf = vdf/sum(vdf(:)); 
tvObj = 0;
for i = 1:numel(neighbors_vdf)
    tvObj = tvObj + sum(sqrt( (vdf(:)-neighbors_vdf{i}(:)).^2 + tvBeta^2 ));
end
f_obj = f_obj + params.gamma*tvObj;
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
    grad = AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,zk),zMask))/bnorm - c;

    % TV regularizer gradient update
    vdf = squeeze(sum(zk,1));
    total = sum(zk(:));
    vdf = vdf/total;
    
    if(sum(vdf(:)) == 0)
        error('VDF is not defined (all zeors)')
    end
    
    if numel(neighbors_vdf) == 2
    	f_tv = [neighbors_vdf{1}; vdf; neighbors_vdf{2}]; 
    elseif imageNum == 1
        f_tv = [vdf; neighbors_vdf{1}]; 
    elseif imageNum == numIms
        f_tv = [vdf; neighbors_vdf{1}];
    end
    
    gradJ = gradientTV(f_tv,tvBeta,D);
    
    for j = 1:numel(vdf)
        gradTV = 0;
        for k = 1:numel(vdf)
            if j == k
                gradTV = gradTV + gradJ(k)*( 1./total - vdf(k)./total );
            else
                gradTV = gradTV - gradJ(k)*vdf(k)./total;
            end     
        end
        grad(:,j) = grad(:,j) + params.gamma*gradTV;
    end
        
    % Backtracking
    stop_backtrack = 0 ;
    while ~stop_backtrack 
        gk = zk - (1/L)*grad ;
        xk = soft(gk,lambda/L) ;
        if isNonnegative
            xk(xk<0) = 0;
        end
        
        % Compute objective at xk
        fit = forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask);
        vdf_xk = squeeze(sum(xk,1));
        vdf_xk = vdf_xk/sum(vdf_xk(:)); 
        tvObj_xk = 0;
        for i = 1:numel(neighbors_vdf)
            tvObj_xk = tvObj_xk + sum(sqrt( (vdf(:)-neighbors_vdf{i}(:)).^2 + tvBeta^2 ));
        end
        
        temp1 = 0.5/bnorm*norm(b(:)-fit(:))^2 + params.gamma*tvObj_xk;
        
        % Compute quadratic approximation at zk
        fit2 = forceMaskToZero(Ax_ft_1D(A0ft_stack,zk),zMask);
        vdf_zk = squeeze(sum(zk,1));
        vdf_zk = vdf_zk/sum(vdf_zk(:)); 
        tvObj_zk = 0;
        for i = 1:numel(neighbors_vdf)
            tvObj_zk = tvObj_zk + sum(sqrt( (vdf(:)-neighbors_vdf{i}(:)).^2 + tvBeta^2 ));
        end
        
        temp2 = 0.5/bnorm*norm(b(:)-fit2(:))^2 +...
            params.gamma*tvObj_zk +...
            (xk(:)-zk(:))'*grad(:) +...
            (L/2)*norm(xk(:)-zk(:))^2;
        
        % Stop backtrack if objective <= quadratic approximation
        if params.noBacktrack
            stop_backtrack = 1 ;
        elseif temp1 <= temp2
            stop_backtrack = 1;
%             params.noBacktrack = 1;
        else
            L = L*beta ;
            if L > 1e15
                keep_going = 0;
                stop_backtrack = 1;
            end
        end
    end
    
    if ~keep_going
        x_hat = xk;
        err = err(1:nIter) ;
        obj = obj(1:nIter) ;
        l_0 = l_0(1:nIter) ;
        return
    end
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k));
    zk = xk + ((t_k-1)/t_kp1)*(xk-xkm1);  

    % Track and display error, objective, sparsity
    prev_f = f_obj;
    f_data = 0.5/bnorm*norm(b-fit)^2;
    f_sparse = lambda * norm(xk(:),1);
    f_tv = params.gamma*tvObj_xk;
    f_obj = f_data + f_sparse + f_tv;   
    err(nIter) = norm(b(:)-fit(:));
    obj(nIter) = f_obj;
    obj1(nIter) = f_data;
    obj2(nIter) = f_sparse;
    obj3(nIter) = f_tv;
    l_0(nIter) = sum(abs(xk(:))>0);
    disp(['Iter ',     num2str(nIter),...
          ' Obj ',     num2str(obj(nIter)),...
          ' L ',       num2str(L),...
          ' ||x||_0 ', num2str(l_0(nIter)),...
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
        plot(fit2);
        title('zk')
        
        subplot(2,3,4)
        fit_gk = forceMaskToZero(Ax_ft_1D(A0ft_stack,gk),zMask);
        plot(fit_gk);
        title('gk')
        
        subplot(2,3,5)
        plot(abs(b-fit));
        title('diff xk')
        
        subplot(2,3,6)
        plot(abs(b-fit2));
        title('diff zk')
        
        pause(0.05)
    end
      
    % Check stopping criterion
    switch stoppingCriterion
        case STOPPING_SUBGRADIENT
            sk = L*(xk-xkm1) +...
                 AtR_ft_2D(A0ft_stack,Ax_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,xk-xkm1),zPad)))/bnorm;
            keep_going = norm(sk(:)) > tolerance*L*max(1,norm(xk(:)));
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
    
    % Update indices
    t_k = t_kp1;
    xkm1 = xk;
end

x_hat = xk;
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