function [x_hat, err, obj] = convADMM_LASSO_TVx_1D(A0ft_stack,b,x_init,x_neighbors,params)
%convADMM_LASSO_1D Image regression by solving LASSO problem 
%                argmin_x 0.5*||Ax-b||^2 + lambda||x||_1
%
% Inputs:
% b          - (n) polar ring image
% A0ft_stack - (n x t) fft of unshifted gaussian basis matrices
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
lambda1 = params.lambda1;
rho1 = params.rho1;
lambda2 = params.lambda2;
rho2 = params.rho2;
maxIter = params.maxIter;
maxGradIters = params.maxGradIters;
isNonnegative = params.isNonnegative;

numIms = params.numIms;
imageNum = params.imageNum;

zPad = params.zeroPad;
zMask = params.zeroMask;
[n,t] = size(A0ft_stack) ;

if ~all(size(x_init)==[n,t])
    error('The dimension of the initial x_init does not match.');
end

b = zeroPad(b,zPad);
bnorm = norm(b(:));
neighbors = numel(x_neighbors);
if neighbors == 2
    D =  [-1  1 0;... 
           0 -1 1];
elseif imageNum == 1
    D = [-1 1];
elseif imageNum == numIms
    D = [1 -1];
end

% Track error and objective
err = nan(1,maxIter);
l1_norm = nan(1,maxIter);
tv_penalty = nan(1,maxIter);
obj = nan(1,maxIter);

% Used to compute gradient
c = AtR_ft_1D(A0ft_stack,b)/bnorm;

% Initialize variables
x_init = forceMaskToZeroArray(x_init,zMask);
xk = x_init;
xkp1 = x_init;

% Arrange neighbors    
if neighbors == 2
    Xk = [x_neighbors{1}(:), xk(:),  x_neighbors{2}(:)]'; 
elseif imageNum == 1
    Xk = [xk(:), x_neighbors{1}(:)]'; 
elseif imageNum == numIms
    Xk = [xk(:), x_neighbors{1}(:)]';
end
Xind = numel(x_neighbors);

yk = xk;
vk = zeros(size(xk));
zk = D*Xk;
uk = zeros(size(zk));

% Initial objective
err(1) = 0.5/bnorm*sum((b-Ax_ft_1D(A0ft_stack,xk)).^2);
l1_norm(1) = lambda1*sum(abs(yk(:)));
tv_penalty(1) = lambda2*sum(abs(zk(:)));
obj(1) = err(1) + l1_norm(1) + tv_penalty(1);

keep_going = 1;
nIter = 1;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;   

    % x-update
    [xkp1,gradIters,params,...
     errX,l1_normX,tv_penaltyX] = convGradDescentTVx(A0ft_stack,xk,b,...
                                        yk,vk,D,Xk,Xind,zk,uk,params );
    % nonegativity
    if isNonnegative
        xkp1(xkp1<0) = 0;
    end
    Xk(Xind,:) = xkp1(:); 
    
    % y-update
    yk = soft(xkp1 + vk/rho1,lambda1/rho1);
    % z-update
    zk = soft(D*Xk + uk/rho2,lambda2/rho2);
    % v-update
    vk = vk + rho1*(xkp1-yk);
    % u-update
    uk = uk + rho2*(D*Xk-zk); 
    
    % objective
    err(nIter) = errX;
    l1_norm(nIter) = l1_normX;
    tv_penalty(nIter) = tv_penaltyX;
    
    obj(nIter) = err(nIter) + l1_norm(nIter) + tv_penalty(nIter);

    disp(['Iter ',      num2str(nIter),...
          ' L ',        num2str(L),...
          ' gradIters ',num2str(gradIters),...
          ' Obj ',      num2str(obj(nIter)),...
          ' RelErr ',   num2str(err(nIter)),...
          ' TV ',       num2str(tv_penalty(nIter)),...
          ' ||x||_1 ',  num2str(l1_norm(nIter)),...
          ' ||x||_0 ',  num2str(sum(xkp1(:) >0))
           ]);
    
    if params.plotProgress
        fit_xkp1 = Ax_ft_1D(A0ft_stack,xkp1);
        figure(1)    
        hold off
        plot(b)
        hold on
        plot(fit_xkp1)
        legend('data','fit')
        
        pause
    end
    
    % Check stopping criterion
    switch stoppingCriterion
        case GRADIENT_NORM
            keep_going = norm(grad(:)) > tolerance;
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
    
    % Update indices
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
