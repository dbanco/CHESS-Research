function X_aligned = apply_perm_to_X(X, K, J, best_perm, shifts)
    % X: size [1, N, J*K]
    % J: number of atoms per dictionary element
    % best_perm: permutation of 1:K
    % shifts: vector of circular shifts of length K
    %
    % Returns X_aligned: permuted and circularly shifted

    [~,N,KJ,T] = size(X);

    X_aligned = zeros(size(X), 'like', X);

    for k = 1:K
        k_idcs = 1+J*(k-1):(J+J*(k-1));  % source indices in X

        X_group = X(:,:,k_idcs,:);    % [1, N, J, T] group
        X_shifted = circshift(X_group, [0, -shift-1, 0, 0]);

        % Store in aligned array
        X_aligned(:,:,tgt_idx,:) = X_shifted;
    end
end

