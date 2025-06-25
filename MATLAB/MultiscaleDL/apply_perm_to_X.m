function X_aligned = apply_perm_to_X(X, J, best_perm, shifts)
    % X: size [1, N, J*K]
    % J: number of atoms per dictionary element
    % best_perm: permutation of 1:K
    % shifts: vector of circular shifts of length K
    %
    % Returns X_aligned: permuted and circularly shifted

    [~, N, JK,T] = size(X);
    K = length(best_perm);
    assert(J*K == JK, 'Third dimension of X must be J*K');
    assert(length(shifts) == K, 'Length of shifts must equal K');

    X_aligned = zeros(size(X), 'like', X);

    for k = 1:K
        src_idx = (best_perm(k)-1)*J + (1:J);  % source indices in X
        tgt_idx = (k-1)*J + (1:J);            % target indices in X_aligned
        X_group = X(:,:,src_idx,:);             % [1, N, J] group
        
        % Store in aligned array
        X_aligned(:,:,tgt_idx,:) = X_group;
    end
end

