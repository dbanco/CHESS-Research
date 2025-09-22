function [Xf, cgst, opt] = x_update_filter_reg_gpu(Xf,Df,bf,opt,N,M,K,J,T)

    if opt.lambda2 > 0
        if opt.ism_init
            Xf = solvedbi_sm(Df, opt.rho1, bf);
        end
        [Xf, cgst, opt] = pcg_filter_reg(Xf,Df,bf,opt,N,M,K,J,T);
    else
        Xf = solvedbi_sm(Df, opt.rho1, bf);
        
        cgst.pit = 1;
    end
    

end
