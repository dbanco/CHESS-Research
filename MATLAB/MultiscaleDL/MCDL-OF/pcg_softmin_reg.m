function [Xf, cgst] = pcg_softmin_reg(Xf, Df, bf, Ddiag, opt, N, M, K, J, T)
% PCG solver for (A^T A + rho I + lambda2 D) X = bf
% with D = diag(Ddiag) in spatial domain (MM surrogate).
%
% Inputs:
%   Xf     : initial iterate in Fourier domain
%   Df     : dictionary in Fourier domain (N+M-1 x ? x K*J)
%   bf     : right-hand side in Fourier domain
%   Ddiag  : spatial diagonal majorizer (same size as X spatial domain)
%   opt    : options struct with fields lambda2, rho1, CGTolX, MaxCGIterX, useGpu
%   N,M,K,J,T : dimensions
%
% Outputs:
%   Xf     : solution in Fourier domain
%   cgst   : struct with PCG status info

lambda2 = opt.lambda2;
rho1    = opt.rho1;
tol     = opt.CGTolX;
mit     = opt.MaxCGIterX;
useGpu  = opt.useGpu;

isn = Xf(:);
xsz = [1, N+M-1, K*J, T];

if useGpu
    isn   = gpuArray(complex(isn));
    Df    = gpuArray(complex(Df));
    Dfh   = gpuArray(complex(conj(Df)));
    bf    = gpuArray(complex(bf));
    Ddiag = gpuArray(Ddiag);
    Aop   = @(uf) sum(pagefun(@times, Df, uf), 3);
    Ahop  = @(uf) pagefun(@times, Dfh, uf);
else
    Dfh   = conj(Df);
    Aop   = @(uf) sum(bsxfun(@times, Df, uf), 3);
    Ahop  = @(uf) bsxfun(@times, Dfh, uf);
end

% A^T A operator in Fourier domain
AhAvop = @(uf) vec(Ahop(Aop(reshape(uf, xsz))));

% Regularizer operator: multiply by D in spatial domain
RegOp = @(uf) vec(fft2(Ddiag .* ifft2(reshape(uf, xsz))));

wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');

% Linear operator for PCG
if lambda2 == 0
    linop = @(uf) AhAvop(uf) + rho1*uf;
else
    linop = @(uf) AhAvop(uf) + rho1*uf + lambda2*RegOp(uf);
end

[Xfvec, flg, rlr, pit, resvec] = pcg(linop, bf(:), tol, mit, [], [], isn);
% disp(resvec);
% semilogy(resvec); xlabel('Iteration'); ylabel('Residual norm');
%     warning(wrn.state, 'MATLAB:ignoreImagPart');

cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);
Xf   = reshape(Xfvec, xsz);

end
