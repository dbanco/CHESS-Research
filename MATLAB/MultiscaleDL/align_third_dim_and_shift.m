function [D_aligned, best_perm, shifts, min_error] = align_third_dim_and_shift(D, Dtrue)
    % D and Dtrue must be of same size, e.g., [1, 45, N] with N <= 5
    assert(isequal(size(D), size(Dtrue)), 'D and Dtrue must have the same size');
    N = size(D, 3);
    
    perms_list = perms(1:N);     % All permutations of 1:N
    num_perms = size(perms_list, 1);
    
    min_error = inf;
    best_perm = [];
    D_aligned = D;
    
    for i = 1:num_perms
        perm = perms_list(i,:);
        D_perm = D(:,:,perm);
        [err, shifts, shifted_D] = compute_shifted_dict_error(D_perm, Dtrue);
        
        if err < min_error
            min_error = err;
            best_perm = perm;
            D_aligned = shifted_D;
        end
    end
    
    fprintf('Best permutation of 3rd dimension: [%s]\n', num2str(best_perm));
    fprintf('Minimum error: %.6f\n', min_error);
end