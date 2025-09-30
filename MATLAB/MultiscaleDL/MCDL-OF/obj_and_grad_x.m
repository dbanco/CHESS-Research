% objective and gradient function (vector -> scalar, vector)
function [f, gvec] = obj_and_grad_x(xvec, Df, Sf, Y, U, rho1, lambda2, K, J)
    % xvec shape -> spatial X
    X = reshape(xvec, size(Y));            % Y is [1,N,KJ,T] spatial
    Xf = fft2(X);
    % forward residual in freq:
    AfXf = sum(bsxfun(@times, Df, Xf), 3); % Aop
    data_term = 0.5 * norm(Sf(:) - AfXf(:))^2;

    % quad term
    quad_term = 0.5 * rho1 * norm(X(:) - Y(:) + U(:))^2;

    % softmin term (you already have compute_softmin)
    R = 0;
    g_soft = zeros(size(X));
    if lambda2 > 0
        [R, g_soft] = compute_softmin(X, K, J);
    end
    f = data_term + quad_term + lambda2 * R;

    % gradient: Atop(AfXf - Sf) + rho*(X - Y + U) + lambda2 * g_soft
    rf = AfXf - Sf;
    grad_data = ifft2(bsxfun(@times, conj(Df), rf),'symmetric'); % returns same shape as X
    g = grad_data + rho1 * (X - Y + U) + lambda2 * g_soft;
    gvec = g(:);
end