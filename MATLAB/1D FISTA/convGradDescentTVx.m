function [xkp1,gradIters,params,err,l1_norm,tv_penalty] = convGradDescentTVx( A0ft_stack,xk,b,yk,vk,D,Xk,Xind,zk,uk,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L = params.L;
beta = params.beta;
rho1 = params.rho1;
rho2 = params.rho2;
maxGradIters = params.maxGradIters;
gradTolerance = params.gradTolerance;
isNonnegative = params.isNonnegative;
zPad = params.zeroPad;
zMask = params.zeroMask;

bnorm = norm(b);
c = AtR_ft_1D(A0ft_stack,b)/bnorm;

for j = 1:maxGradIters
    % Compute gradient of f
    grad = AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask))/bnorm - c +...
           rho1*(xk - yk + vk) + ...
           rho2*reshape(D(:,Xind)'*(D*Xk - zk + uk),size(xk));
    
       % backtracking step size search
    stop_backtrack = 0;
    while ~stop_backtrack 
        xkp1 = xk - (1/L)*grad;
        Xkp1 = Xk;
        Xkp1(Xind,:) = xkp1(:); 
        % Compute objective at xk/xkp1 and criterion to stop backtracking
        fit_xk = forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask);
        f_xk = 0.5/bnorm*sum((b-fit_xk).^2) +...
            rho1/2*sum(( xk(:) - yk(:) + vk(:) ).^2)+...
            rho2/2*sum(sum(sum(( D*Xk - zk + uk).^2)));
        fit_xkp1 = forceMaskToZero(Ax_ft_1D(A0ft_stack,xkp1),zMask);
        
        err = 0.5/bnorm*sum((b-fit_xkp1).^2);
        l1_norm = rho1/2*sum(( xkp1(:) - yk(:) + vk(:) ).^2);
        tv_penalty = rho2/2*sum(sum(sum(( D*Xkp1 - zk + uk ).^2)));   
        f_xkp1 = err + l1_norm + tv_penalty;     
        
        criterion = f_xk - (0.3/L)*sum(grad(:).^2);

        % Stop backtrack logic
        if params.noBacktrack
            stop_backtrack = 1 ;
        elseif f_xkp1 <= criterion
            stop_backtrack = 1;
            params.noBacktrack = 1;
        else
            L = L*beta ;
            if L > 1e50
                keep_going = 0;
                stop_backtrack = 1;
            end
        end
    end
    if norm(grad(:)) < gradTolerance
        break
    end
    xk = xkp1;
    Xk(Xind,:) = xk(:);
end
gradIters = j;
params.L = L;
end

