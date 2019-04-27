function [x_hat, err, t_k, L, obj, l_0] = space_wasserstein_FISTA_Circulant(A0ft_stack,b,neighbors_vdf,D,x_init,params)
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
wLam = params.wLam;
beta = params.beta;
maxIter = params.maxIterReg;
isNonnegative = params.isNonnegative;
zPad = params.zeroPad;
zMask = params.zeroMask;
[m,n,t,r] = size(A0ft_stack) ;

if ~all(size(x_init)==[m,n,t,r])
    error('The dimension of the initial xk does not match.');
end

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
obj1 = nan(1,maxIter);
obj2 = nan(1,maxIter);
obj3 = nan(1,maxIter);
l_0 = nan(1,maxIter);

% Initial sparsity and objective
f = 0.5*norm(b-Ax_ft_2D(A0ft_stack,x_init))^2 +...
    lambda * norm(x_init(:),1);

% Add entropic reg wasserstein distance vdf term  
vdf = squeeze(sum(sum(x_init,1),2));
vdf = vdf/sum(vdf(:)); 
wObj = WassersteinObjective(vdf(:), neighbors_vdf(:), wLam, D);
f = f + 0.5*params.gamma*wObj;

% Used to compute gradient
c = AtR_ft_2D(A0ft_stack,b);

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
    grad = AtR_ft_2D(A0ft_stack,forceMaskToZero(Ax_ft_2D(A0ft_stack,zk),zMask)) - c;

    % Wasserstein regularizer gradient update
    vdf = squeeze(sum(sum(zk,1),2));
    vdf = vdf/sum(vdf(:));
    if(sum(vdf(:)) == 0)
        error('VDF is not defined (all zeors)')
    end
    gradW = zeros(t*r,1);
    wObj_zk = 0;
    for i = 1:numel(neighbors_vdf)
        [gradWi, Wdi] = WassersteinGrad( vdf(:), neighbors_vdf{i}(:), wLam, D );
        gradW = gradW + gradWi;
        wObj_zk = wObj_zk + Wdi;
    end
   gradW = reshape(gradW,[t,r])./numel(neighbors_vdf);
    for i = 1:t
        for j = 1:r
            grad(:,:,i,j) = grad(:,:,i,j) + ...
            params.gamma*(gradW(i,j)./sum(zk(:)) - sum(gradW(:).*vdf(:))./sum(zk(:)));
        end
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
        fit = forceMaskToZero(Ax_ft_2D(A0ft_stack,xk),zMask);
        vdf_xk = squeeze(sum(sum(xk,1),2));
        vdf_xk = vdf_xk/sum(vdf_xk(:)); 
        wObj_xk = WassersteinObjective(vdf_xk(:),neighbors_vdf(:),wLam,D);
        temp1 = 0.5*norm(b(:)-fit(:))^2 + 0.5*params.gamma*wObj_xk;
        
        % Compute quadratic approximation at zk
        fit2 = forceMaskToZero(Ax_ft_2D(A0ft_stack,zk),zMask);
        vdf_zk = squeeze(sum(sum(zk,1),2));
        vdf_zk = vdf_zk/sum(vdf_zk(:)); 
        temp2 = 0.5*norm(b(:)-fit2(:))^2 +...
            0.5*params.gamma*wObj_zk +...
            (xk(:)-zk(:))'*grad(:) +...
            (L/2)*norm(xk(:)-zk(:))^2;
        
        % Stop backtrack if objective <= quadratic approximation
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*beta ;
            if L > 1e50
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
    prev_f = f;
    f_data = 0.5*norm(b-fit)^2;
    f_sparse = lambda * norm(xk(:),1);
    f_wasserstein = 0.5*params.gamma*wObj_xk;
    f = f_data + f_sparse + f_wasserstein;   
    err(nIter) = norm(b(:)-fit(:))/norm(b(:));
    obj(nIter) = f;
    obj1(nIter) = f_data;
    obj2(nIter) = f_sparse;
    obj3(nIter) = f_wasserstein;
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
        imshow(b,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('img')
        
        subplot(2,3,2)
        imshow(Ax_ft_2D(A0ft_stack,xk),'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('xk')
        
        subplot(2,3,3)
        imshow(fit2,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('zk')
        
        subplot(2,3,4)
        fit_gk = forceMaskToZero(Ax_ft_2D(A0ft_stack,gk),zMask);
        imshow(fit_gk,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('gk')
        
        subplot(2,3,5)
        imshow(abs(b-fit),'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('diff xk')
        
        subplot(2,3,6)
        imshow(abs(b-fit2),'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('diff zk')
        
        pause(0.05)
    end
      
    % Check stopping criterion
    switch stoppingCriterion
        case STOPPING_SUBGRADIENT
            sk = L*(xk-xkm1) +...
                 AtR_ft_2D(A0ft_stack,Ax_ft_2D(A0ft_stack,forceMaskToZero(Ax_ft_2D(A0ft_stack,xk-xkm1),zPad)));
            keep_going = norm(sk(:)) > tolerance*L*max(1,norm(xk(:)));
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f-prev_f);
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