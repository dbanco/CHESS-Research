function [x, cgst] = x_update_switch(Xf,D,rho,b,opt,N,M,K,J,T,lambda2)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
switch opt.regularizer
    case 'filter1'
        [x, cgst] = x_update_filter_reg_gpu(Xf,D,rho,b,opt,N,M,K,J,T,lambda2);
    case 'soft-min flat'
        [x, cgst] = x_update_softmin_reg_gpu(Xf,D,rho,b,opt,N,M,K,J,T,lambda2);
end

end