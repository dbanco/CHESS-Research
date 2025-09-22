function reg = compute_regularizer(X,K,J,regularizer,transpose)

if nargin < 5
    transpose = false;
end
switch regularizer
    case 'filter1'
        kernel = [[ 0  0  0;
                    0  1  0;
                    0  0  0],...
                  [-1 -1 -1;
                   -1 -1 -1;
                   -1 -1 -1]];

        reg = compute_filter_regularizer(kernel,X,K,J,transpose);
    case 'soft-min flat'
        [reg,~] = compute_time_reg_softmin_flat(X_flat, K, J, 1e-3);
end

end