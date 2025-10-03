function [Xf, cgst, opt] = x_update_softmin(X0,Df,bf,Sf,Y,U,opt,N,M,K,J,T)

switch opt.optimizer
    case 'LBFGS'
        [Xf, cgst, opt] = x_update_softmin_lbfgs_gpu(X0,Df,bf,Sf,Y,U,opt,N,M,K,J,T);
    case 'LinMM'
        [Xf, cgst, opt] = x_update_softmin_mm_ism(X0,Df,bf,Sf,Y,U,opt,N,M,K,J,T);
    case 'QuadMM'
        [Xf, cgst, opt] = x_update_softmin_mm_pcg_gpu(X0,Df,bf,Sf,Y,U,opt,N,M,K,J,T);
    case 'FISTA'
        [Xf, cgst, opt] = x_update_softmin_fista(X0,Df,bf,Y,U,opt,K,J,T);
end