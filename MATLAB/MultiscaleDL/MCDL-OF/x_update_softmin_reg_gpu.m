function [Xf, cgst, opt] = x_update_softmin_reg_gpu(X,Df,bf,S,Y,U,opt,N,M,K,J,T,Atop)
    opts.b = 0;
    opts.L = opt.L;

    if opt.lambda2 > 0
        if opt.ism_init
            Xf = solvedbi_sm(Df, opt.rho1, bf);
            X = ifft2(Xf,'symmetric');
        end
        [X, objVals,L] = grad_desc_softmin_x_update(X,Df,S,Y,U,opt,Atop,K,J,M,opts);
        Xf = fft2(X);
        cgst.pit = numel(objVals);
        opt.L = L;
    else
        Xf = solvedbi_sm(Df, opt.rho1, bf);
        cgst.pit = 1;
    end

end