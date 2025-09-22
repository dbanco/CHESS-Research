function [Xf, cgst, opt] = pcg_filter_reg(Xf,Df,bf,opt,N,M,K,J,T)

lambda2 = opt.lambda2;
rho1 = opt.rho1;
tol = opt.CGTolX;
mit = opt.MaxCGIterX;
useGpu = opt.useGpu;
isn = Xf(:);

xsz = [1,N+M-1,K*J,T];

if useGpu
    isn = gpuArray(complex(isn));
    Df = gpuArray(complex(Df));
    Dfh = gpuArray(complex(conj(Df)));
    bf = gpuArray(complex(bf));
    Aop = @(uf) sum(pagefun(@times, Df, uf), 3);
    Ahop = @(uf) pagefun(@times, Dfh, uf);
else
    Dfh = conj(Df);
    Aop = @(uf) sum(bsxfun(@times, Df, uf), 3);
    Ahop = @(uf) bsxfun(@times, Dfh, uf);
end

AhAvop = @(uf) vec(Ahop(Aop(reshape(uf, xsz))));
AhAvop2 = @(uf) AtAx_filter(uf,xsz,K,J,opt);

wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');

if (lambda2 == 0)
    [Xfvec,flg,rlr,pit,~] = pcg(@(uf) AhAvop(uf) + rho1*uf,...
    bf(:), tol, mit, [], [], isn);
else
    [Xfvec,flg,rlr,pit,~] = pcg(@(uf) AhAvop(uf) + rho1*uf + lambda2*AhAvop2(uf),...
    bf(:), tol, mit, [], [], isn);
end

warning(wrn.state, 'MATLAB:ignoreImagPart'); % 
cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);

Xf = reshape(Xfvec, xsz);

end

function AtAxf = AtAx_filter(uf,xsz,K,J,opt)

    X = ifft2(reshape(uf,xsz),'symmetric');
    X1 = compute_regularizer(X,K,J,opt.regularizer);
    X = compute_regularizer(X1,K,J,opt.regularizer,true);
    
    AtAxf = vec(fft2(X));

end