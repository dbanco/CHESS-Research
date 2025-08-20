function [D_aligned, d_shifts] = align_third_dim_and_shift(D, Dtrue)
    % D and Dtrue must be of same size, e.g., [1, 45, N] with N <= 5
    [~, N, K] = size(D);
    assert(isequal(size(D), size(Dtrue)), 'D and Dtrue must be the same size');
    
    D_aligned = D;
    d_shifts = zeros(K,1);
    
    for k = 1:K
        d = squeeze(D(1,:,k));
        d_true = squeeze(Dtrue(1,:,k));
        corr = xcorr(d, d_true, 'coeff');
        [~, max_idx] = max(corr);
        shift = mod(max_idx - N - 1, N);   % shift to align d to d_true
        d_shifted = circshift(d, [0, -shift-1]);
        D_aligned(1,:,k) = d_shifted;
        d_shifts(k) = -shift-1;
    end
end