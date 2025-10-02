function [Xf, cgst, opt] = x_update_softmin_mm_pcg_gpu(X0,Df,bf,Sf,Y,U,opt,N,M,K,J,T)

if opt.lambda2 > 0
    if opt.ism_init
        Xf = solvedbi_sm(Df, opt.rho1, bf);
        X0 = ifft2(Xf,'symmetric');
    else
        Xf = fft2(X0);
    end
    
    [gradR, Ddiag] = softmin_grad_majorizer(X0, K, J, opt.tau);
    Ddiag(Ddiag < 1e-6) = 1e-6;
    bf_eff = bf - opt.lambda2*fft2(-gradR + Ddiag.*X0);
    [Xf, cgst] = pcg_softmin_reg(Xf, Df, bf_eff, Ddiag, opt, N, M, K, J, T);
else
    Xf = solvedbi_sm(Df, opt.rho1, bf);
    cgst.pit = 1;
end
    

end