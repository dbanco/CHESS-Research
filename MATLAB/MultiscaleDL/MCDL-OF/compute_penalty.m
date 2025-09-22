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
        Jreg = compute_filter_regularizer_penalty(kernel,X,K,J);
    case 'filter3'
        k_spatial = fspecial('gaussian', [3 3], 1);
        kernel = zeros(3,3,2);
        kernel(:,:,1) = -k_spatial;
        kernel(:,:,2) = k_spatial;
        Jreg = compute_filter_regularizer_penalty(kernel,X,K,J);
    case 'softmin'
        Jreg = compute_softmin(X,K,J); 
end

end