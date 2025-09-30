Xtest = fft2(randn(size(X)));  % same shape as X
Rtest = fft2(randn(size(Sf)));  % same shape as S

% Apply forward and adjoint
AX = Aop(Xtest);          % forward
AtR = Atop(Rtest);              % adjoint

Aop = @(xf) sum(bsxfun(@times,Df,xf),3);
Atop = @(rf) bsxfun(@times, conj(Df), rf);

% Compute inner products
lhs = sum(conj(AX(:)) .* Rtest(:));  % <A x, y>
rhs = sum(conj(Xtest(:)) .* AtR(:)); % <x, A* y>

fprintf('Adjoint test: |lhs-rhs| / |lhs| = %.3e\n', abs(lhs-rhs)/abs(lhs));