function [Xf, cgst, opt] = x_update_switch(X,Xf,Df,bf,S,Y,U,opt,N,M,K,J,T,Atop)
%x_update_switch  Update sparse codes with different regularization options.
%
%   [Xf, cgst, opt] = x_update_switch(X,Xf,Df,bf,S,Y,U,opt,N,M,K,J,T,Atop)
%
%   This function selects and runs the appropriate X-update routine in an
%   ADMM or dictionary learning loop, depending on the regularizer specified
%   in OPT.regularizer. It provides a unified interface for different
%   regularization strategies.
%
%   INPUTS:
%     X       - Current estimate of the coefficient maps in the spatial domain.
%     Xf      - Current estimate of the coefficient maps in the frequency domain.
%     Df      - Dictionary filters in the frequency domain.
%     bf      - Right-hand side term in the normal equations (frequency domain).
%     S       - Observed signal/data (spatial domain).
%     Y, U    - ADMM auxiliary and dual variables, respectively.
%     opt     - Options struct with fields:
%                  .regularizer : String specifying the regularization type.
%                                 Supported values:
%                                   'filter1', 'filter2', 'filter3'
%                                      → L2 filter-based regularization
%                                   'softmin'
%                                      → Log-sum-exp softmin regularization
%     N, M    - Problem dimensions for coefficient maps.
%     K, J    - Number of filters/scales or related indexing parameters.
%     T       - Number of time frames (for spatiotemporal problems).
%     Atop    - Function handle for adjoint operator A^T (used in softmin case).
%
%   OUTPUTS:
%     Xf      - Updated coefficient maps in the frequency domain.
%     cgst    - Struct with conjugate gradient solver status:
%                  .flg  : Convergence flag from PCG
%                  .rlr  : Relative residual norm
%                  .pit  : Number of iterations performed
%     opt     - Updated options struct (may contain updated tolerances, etc.).
%
%   NOTES:
%     • 'filter1', 'filter2', 'filter3' currently call the same update routine
%       (x_update_filter_reg_gpu) but may correspond to different kernel choices.
%     • 'softmin' uses x_update_softmin_reg_gpu, which operates in the spatial
%       domain and applies a log-sum-exp type regularizer.
%     • Xf is always returned in the frequency domain for consistency.
%
%   See also: x_update_filter_reg_gpu, x_update_softmin_reg_gpu

switch opt.regularizer
    case 'filter1'
        [Xf, cgst, opt] = x_update_filter_reg_gpu(Xf,Df,bf,opt,N,M,K,J,T);
    case 'filter2'
        [Xf, cgst, opt] = x_update_filter_reg_gpu(Xf,Df,bf,opt,N,M,K,J,T);
    case 'filter3'
        [Xf, cgst, opt] = x_update_filter_reg_gpu(Xf,Df,bf,opt,N,M,K,J,T);
    case 'softmin'
        [Xf, cgst, opt] = x_update_softmin_reg_gpu(X,Df,bf,S,Y,U,opt,N,M,K,J,T,Atop);
        
end

end