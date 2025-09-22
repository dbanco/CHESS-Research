function reg = compute_regularizer(X,K,J,regularizer,transpose)

if nargin < 5
    transpose = false;
end
switch regularizer
    case 'filter1'
        kernel = zeros(3,3,2);
        kernel(:,:,1) = [ 0  0  0;
                          0  1  0;
                          0  0  0];
        kernel(:,:,2) = [-1 -1 -1;
                         -1 -1 -1;
                         -1 -1 -1];
        reg = compute_filter_regularizer(kernel,X,K,J,transpose);
    case 'filter2'
        kernel = zeros(3,3,3);
        kernel(:,:,1) = 0.5*[-1 -1 -1;
                             -1 -1 -1;
                             -1 -1 -1];
        kernel(:,:,2) = [ 0  0  0;
                          0  1  0;
                          0  0  0];
        kernel(:,:,3) = 0.5*[-1 -1 -1;
                             -1 -1 -1;
                             -1 -1 -1];
        reg = compute_filter_regularizer(kernel,X,K,J,transpose);
    case 'filter3'
        k_spatial = fspecial('gaussian', [3 3], 1);
        kernel = zeros(3,3,2);
        kernel(:,:,1) = -k_spatial;
        kernel(:,:,2) = k_spatial;
        reg = compute_filter_regularizer(kernel,X,K,J,transpose);
    case 'soft-min flat'
        [reg,~] = compute_softmin(X_flat, K, J, 1e-3);
end

end