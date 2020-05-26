function [xkp1,gradIters,params,err,l1_norm] = convGradDescent( A0ft_stack,xk,b,yk,vk,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L = params.L;
beta = params.beta;
rho = params.rho;
maxGradIters = params.maxGradIters;
gradTolerance = params.gradTolerance;
isNonnegative = params.isNonnegative;
zPad = params.zeroPad;
zMask = params.zeroMask;

bnorm = norm(b);
c = AtR_ft_1D(A0ft_stack,b)/bnorm;

for j = 1:maxGradIters
    grad = AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask))/bnorm - c +...
           rho*(xk - yk + vk);

    % backtracking step size search
    stop_backtrack = 0;
    while ~stop_backtrack 
        xkp1 = xk - (1/L)*grad;

        % Compute objective at xk/xkp1 and criterion to stop backtracking
        fit_xk = forceMaskToZero(Ax_ft_1D(A0ft_stack,xk),zMask);
        f_xk = 0.5/bnorm*sum((b-fit_xk).^2) + rho/2*sum(( xk(:) - yk(:) + vk(:) ).^2);
        fit_xkp1 = forceMaskToZero(Ax_ft_1D(A0ft_stack,xkp1),zMask);
        err = 0.5/bnorm*sum((b-fit_xkp1).^2);
        l1_norm = rho/2*sum(( xkp1(:) - yk(:) + vk(:) ).^2);
        f_xkp1 = err + l1_norm;
        
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
end
gradIters = j;
params.L = L;
end

