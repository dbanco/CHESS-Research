function [total_error, shifts, aligned_D] = compute_shifted_dict_error(D, Dtrue)
    % Both D and Dtrue are [1, N, K]
    % Returns the total error after best circular alignment per atom
    % along with the shift indices used and aligned D

    [~, N, K] = size(D);
    assert(isequal(size(D), size(Dtrue)), 'D and Dtrue must be the same size');

    total_error = 0;
    shifts = zeros(1, K);
    aligned_D = zeros(size(D), 'like', D);

    for k = 1:K
        d = squeeze(D(1,:,k));
        d_true = squeeze(Dtrue(1,:,k));
        
        % Cross-correlation to find shift (assumes circular shift)
        corr = xcorr(d, d_true, 'coeff');  % normalized cross-correlation
        [~, max_idx] = max(corr);
        shift = mod(max_idx - N - 1, N);   % shift to align d to d_true

        % Apply circular shift
        d_shifted = circshift(d, [0, -shift]);
        aligned_D(1,:,k) = d_shifted;

        % Compute error
        err = norm(d_shifted - d_true)^2;
        total_error = total_error + err;
        shifts(k) = shift;
    end

    total_error = sqrt(total_error);  % Frobenius-like total error
end
