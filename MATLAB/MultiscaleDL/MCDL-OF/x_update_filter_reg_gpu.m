function [x, cgst] = x_update_filter_reg_gpu(Xf,D,rho,b,opt,N,M,K,J,T,lambda2)

tol = opt.CGTolX;
mit = opt.MaxCGIter;
isn = Xf;
useGpu = opt.useGpu;

xsz = [1,N+M-1,K*J,T];

if useGpu
    D = gpuArray(complex(D));
    ah = gpuArray(complex(conj(D)));
    b = gpuArray(complex(b));
    Aop = @(u) ifft2(sum(pagefun(@times, D, u), 3),'symmetric');
    Ahop = @(u) pagefun(@times, ah, u);
else
    ah = conj(D);
    Aop = @(u) ifft2(sum(bsxfun(@times, D, u), 3),'symmetric');
    Ahop = @(u) bsxfun(@times, ah, u);
end

AhAvop = @(u) vec(Ahop(fft2(Aop(reshape(u, xsz))))); % This is where maskPad was
AhAvop2 = @(u) AtAx_filter(u,xsz,K,J,opt);

wrn = warning('query','MATLAB:ignoreImagPart');
warning('off', 'MATLAB:ignoreImagPart');

if (lambda2 == 0)
    [xv,flg,rlr,pit,resvec] = pcg(@(u) AhAvop(u) + rho*u,...
    b(:), tol, mit, [], [], isn);
else
    [xv,flg,rlr,pit,resvec] = pcg(@(u) AhAvop(u) + rho*u + lambda2*AhAvop2(u),...
    b(:), tol, mit, [], [], isn);
end

warning(wrn.state, 'MATLAB:ignoreImagPart'); % 
cgst = struct('flg', flg, 'rlr', rlr, 'pit', pit);

x = reshape(xv, xsz);

end

function AtAxf = AtAx_filter(u,xsz,K,J,opt)

X = real(ifft2(reshape(u,xsz),'symmetric'));

X1 = compute_regularizer(X,K,J,opt.regularizer);
X = compute_regularizer(X1,K,J,opt.regularizer,true);

AtAxf = vec(fft2(X));



end
