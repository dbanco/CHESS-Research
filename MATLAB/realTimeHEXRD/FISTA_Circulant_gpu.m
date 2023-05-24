function [Xk, err, obj, l_1, t_k, L] = FISTA_Circulant_gpu(Df,b,params)
%FISTA_Circulant Image regression by solving LASSO problem 
%                argmin_x 0.5*||Ax-b||^2 + lambda||x||_1
%
%   Implementation of Fast Iterative Shrinkage-Thresholding Algorithm using 
%   convolutional subroutines for circulant matrix computations as
%   described in:
%   A. Beck and M. Teboulle, â€œA Fast Iterative Shrinkage-Thresholding 
%       Algorithm for Linear Inverse Problems,â€? SIAM Journal on Imaging 
%       Sciences, vol. 2, no. 1, pp. 183202, Jan. 2009.
%
% Inputs:
% b          - (N x M) polar ring image
% Df - (N x M x K x T) fft2 of unshifted gaussian basis matrices
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
% x_hat - (N x M x K x T) solution 
% err - (nIters) relative error of solution at each iteration
% obj - (nIters) objective value of solution at each iteration
% l_0 - (nIters) sparsity of solution at each iteration


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

% [N,M,K,T] = size(Df);
% if ~all(size(x_init)==[N,M,K,T])
%     error('The dimension of the initial xk does not match.');
% end

Df = gpuArray(complex(Df));
bf = gpuArray(complex(fft2(b)));
c = ifft2(Atx_gpu(Df,bf),'symmetric');
Dtbf = Atx_gpu(Df,bf);

% Initialize solution
% Xkf = solvedbi_sm(Df, 1, Dtbf);
% Xk = ifft2(Xkf, 'symmetric');
Xk = gpuArray(zeros(size(Dtbf)));
Xkf = gpuArray(complex(zeros(size(Dtbf))));
Xkf = Xk;
Xkm1 = Xk;
Zk = Xk;
Zkf = Xkf;

% Track error, objective, and sparsity
err = nan(1,maxIter);
obj = nan(1,maxIter);
l_1 = nan(1,maxIter);
l_0 = nan(1,maxIter);
if params.verbose
    disp('Iter    Obj     L        ||x||_0         ||x||_1     RelErr ')
    disp_string = '%i       %0.2f   %0.1f     %i     %0.2f     %0.2f \n';
end


% Initial sparsity and objective
f = 0.5*norm(b-ifft2(Ax_gpu(Df,Xkf),'symmetric'))^2 +...
    lambda * norm(Xk(:),1);

% Used to compute gradient
t_k = 1;
keep_going = 1;
nIter = 0;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;        
        
    % Compute gradient of f
    grad = ifft2(Atx_gpu(Df,Ax_gpu(Df,Zkf)),'symmetric') - c; % gradient of f at zk
    
    % Backtracking Line Search
    stop_backtrack = 0 ;
    while ~stop_backtrack 
        
        %l1/nonnegative-proximal
        gk = Zk - (1/L)*grad ;
        Xk = soft(gk,lambda/L) ;
        if isNonnegative
            Xk(Xk<0) = 0;
        end
        Xkf = fft2(Xk);
        
        % Compute quadratic approximation at yk
        fitX = ifft2(Ax_gpu(Df,Xkf),'symmetric');
        fitZ = ifft2(Ax_gpu(Df,Zkf),'symmetric');
        
        temp1 = 0.5*sum((b(:)-fitX(:)).^2)  + lambda*sum(abs(Xk(:)));
        temp2 = 0.5*sum((b(:)-fitZ(:)).^2) + lambda*sum(abs(Xk(:))) +...
            (Xk(:)-Zk(:))'*grad(:) + (L/2)*norm(Xk(:)-Zk(:))^2;
        
        % Stop backtrack if objective <= quadratic approximation
        if temp1 <= temp2
            stop_backtrack = 1 ;
        elseif params.noBacktrack
            stop_backtrack = 1;
        else
            L = L*beta ;
        end
    end
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k));
    Zk = Xk + ((t_k-1)/t_kp1)*(Xk-Xkm1);    
    Zkf = fft2(Zk);

    % Track and display error, objective, sparsity
    prev_f = f;
    err(nIter) = norm(b(:)-fitX(:))^2;
    l_1(nIter) = sum(abs(Xk(:)));
    l_0(nIter) = sum(abs(Xk(:))>eps*10);
    f = 0.5*err(nIter) + lambda*l_1(nIter);
    err(nIter) = norm(b(:)-fitX(:))/norm(b(:));
    obj(nIter) = f;

    if params.verbose
        fprintf(disp_string,nIter,f,L,l_0(nIter),l_1(nIter),err(nIter))
    end
    
    if params.plotProgress
        lim1 = 0;
        lim2 = max(b(:));
        figure(1)
       
        subplot(2,3,1)
        imshow(b,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('img')
        
        subplot(2,3,2)
        imshow(ifft2(Ax_gpu(Df,Xkf),'symmetric'),'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('xk')
        
        subplot(2,3,3)
        imshow(fitZ,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('zk')
        
        subplot(2,3,4)
        fit_gk = ifft2(Ax_gpu(Df,fft2(gk)),'symmetric');
        imshow(fit_gk,'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('gk')
        
        subplot(2,3,5)
        imshow(abs(b-fitX),'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('diff xk')
        
        subplot(2,3,6)
        imshow(abs(b-fitZ),'DisplayRange',[lim1 lim2],'Colormap',jet);
        title('diff zk')
        
        pause(0.05)
    end
    

    % Check stopping criterion
    switch stoppingCriterion
        case 'SUBGRADIENT'
            sk = L*(Xk-Xkm1) +...
                 Atx_gpu(Df,Ax_gpu(Df,fft2(Xk-Xkm1)));
            keep_going = norm(sk(:)) > tolerance*L*max(1,norm(Xk(:)));
        case 'OBJECTIVE_VALUE'
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            criterionObjective = abs(f-prev_f)/abs(f);
            keep_going =  (criterionObjective > tolerance);
        case 'COEF_CHANGE'
            diff_x = sum(abs(Xk(:)-Xkm1(:)))/numel(Xk);
            keep_going = (diff_x > tolerance);
        otherwise
            error('Undefined stopping criterion.');
    end
    
    % Update indices
    t_k = t_kp1;
    Xkm1 = Xk;
end

err = err(1:nIter) ;
obj = obj(1:nIter) ;
l_1 = l_1(1:nIter) ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end