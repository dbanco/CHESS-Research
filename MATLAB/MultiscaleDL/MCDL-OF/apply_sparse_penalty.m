function [Y] = apply_sparse_penalty(X,lambda,k,opt,a_n)
        
    switch opt.Penalty
        case 'l1-norm'
            Y = shrink(X, (lambda/opt.rho1)*opt.L1Weight);
        case 'log'
            if k <= opt.l1_iters
                Y = shrink(X, (lambda/opt.rho1));
            else
                Y = log_shrink(X, (lambda/opt.rho1),a_n);
            end
    end
    if opt.NonNegCoef
        Y(Y < 0) = 0;
    end
    
    if opt.NoBndryCross
        %Y((end-max(dsz(1,:))+2):end,:,:,:) = 0;
        Y((end-size(D0,1)+2):end,:,:,:) = 0;
        %Y(:,(end-max(dsz(2,:))+2):end,:,:) = 0;
        Y(:,(end-size(D0,2)+2):end,:,:) = 0;
    end

end


function u = shrink(v, lambda)
    if isscalar(lambda)
        u = sign(v).*max(0, abs(v) - lambda);
    else
        u = sign(v).*max(0, bsxfun(@minus, abs(v), lambda));
    end
end

function u = log_shrink(v, lambda, a)
    low_v = v < lambda;
    y1 = abs(v)/2 - 1/(2*a);
    y2 = abs(v)/2 + 1/(2*a);
    disc = y2.^2 - lambda/a;
    disc = max(disc, 0);
    u = sign(v).*(y1 + sqrt(disc));
    u(low_v) = 0;
end
