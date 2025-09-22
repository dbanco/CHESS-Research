% pick small random test arrays in spatial domain
opt.regularizer = 'filter1';

K = 1;
J = 8;
N = 20;
T = 8;
rng(0);
X = randn(1,N, K*J, T);   % adapt sizes to match compute_regularizer shape [K J N T]
Y = randn(size(compute_regularizer(X, K, J, opt.regularizer))); % shape of R*X

% compute forward and adjoint
RX = compute_regularizer(X, K, J, opt.regularizer);
RTY = compute_regularizer(Y, K, J, opt.regularizer,true);

lhs = sum(conj(RX(:)).* Y(:));  % <R x, y>
rhs = sum(conj(X(:)).* RTY(:)); % <x, R^T y>

rel_err = abs(lhs - rhs) / max(1e-12, abs(lhs));
fprintf('Adjoint test (R): rel error = %.3e (lhs=%.3e rhs=%.3e)\n', rel_err, lhs, rhs);
