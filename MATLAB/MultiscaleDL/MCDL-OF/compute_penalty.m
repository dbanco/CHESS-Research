function Jreg = compute_penalty(X,K,J,regularizer)

switch regularizer
    case 'filter1'
        kernel = [[ 0  0  0;
                    0  1  0;
                    0  0  0],...
                  [-1 -1 -1;
                   -1 -1 -1;
                   -1 -1 -1]];
        Jreg = compute_filter_regularizer_penalty(kernel,X,K,J);
    case 'soft-min flat'
        Jreg = compute_time_reg_softmin_flat(X,K,J);
end

end