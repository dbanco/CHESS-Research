function Jl1 = compute_sparse_penalty(Y,a_n,opt,k)

switch opt.Penalty
    case 'l1-norm'
        Jl1 = sum(abs(vec(bsxfun(@times, opt.L1Weight, Y))));
    case 'log'
        if k <= opt.l1_iters
            Jl1 = sum(abs(vec(bsxfun(@times, opt.L1Weight, Y))));
        else
            Jl1 = sum(vec(log(1 + a_n.*abs(Y))/a_n));
        end
end
    
end