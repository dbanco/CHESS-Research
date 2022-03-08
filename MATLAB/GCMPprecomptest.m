function a = GCMPprecomptest(A0ft_stack, x, params)
%GCMP Group Convolutional Matching Pursuit
% Solves a relaxed version of the convolutional sparse coding problem 
%    min ||a||_0,inf
%    s.t. x = Da 
% referred to as the P^epsilon_0,inf problem
%    min ||a||_0,inf
%    s.t. ||x - Da||^2 <= epsilon
% where generally speaking,
% D is a BCCB dictionary (block circulant with circulant blocks)
% a are linear coefficients defing the solution
% x is input data to be coded
% epsilon is an error bound.
%
% Implemented mostly as described in
% A Greedy Approach to l_0, Based Convolutional Sparse Coding
% by Elad Plaut and Raja Giryes
%
%  Different in that we will use 
%  -||x - Da||/||x|| instead of the MSE ||x - Da||^2
%  -Implement overlap detection by circular convolutions
%  -Code makes below assumptions on the dimensions of data/dictionary
%
% Inputs:
% x - (m x n) polar ring image
% A0ft_stack - (m x n x t x r) fft2 of unshifted gaussian basis matrices
% params - struct containing the following field
%   epsilon - relative error stopping parameter in (0,1)
%   delta - change in relative error stopping parameter in (0,1)
%   zeroMask - (#pixels x 2) indices of unobserved pixels 
%   isNonnegative - flag to enforce nonnegative solution
%   showImage - flag to display image as each spot is placed
%
% Outputs:
% a - (m x n x t x r) solution 

a = zeros(size(A0ft_stack));
R = x;
k = 0;
[m,n,t,r] = size(A0ft_stack);

zMask = params.zeroMask;
if params.showImage
    uplim = max(x(:));
    figure(1)
    subplot(2,1,1)
    imshow(x,'DisplayRange',[0 uplim],'Colormap',jet)
    title('x')
    pause(0.05)
end

% Precompute norm of data
x_norm = norm(x(:));

% Compute overlap matrices
[overlapMinRow,overlapMaxRow,overlapMinCol,overlapMaxCol] = overlapMatrix(A0ft_stack,t,r);

while norm(R(:))/x_norm > params.epsilon
    % Compute inner products
    b = AtR_ft_2D(A0ft_stack,R);
    
    % Enforce nonnegativity
    if params.isNonnegative
        b(b<1e-12) = 0;
    end
    
    % Find max inner product
    [bi_star, i_star] = max(b(:));
    
    while bi_star > 1e-12

        % Update solution
        a(i_star) = a(i_star) + bi_star;
        
        if params.showImage
            fit = forceMaskToZero(Ax_ft_2D(A0ft_stack,a),zMask);
            figure(1)
            subplot(2,1,2)
            imshow(fit,'DisplayRange',[0 uplim],'Colormap',jet)
            title('Da')
            pause(0.05)
        end

        % Convert linear index i_star to 4D dictionary index
        [row_star,col_star,t_star,r_star] = ind2sub([m,n,t,r],i_star);
        
        % Loop over all dictionary atom shaps to remove overlapping atoms
        fprintf('[-')
        for i = 1:t
            fprintf('-%i-',round(100*i/t))
            for j = 1:r
                DrowMin = overlapMinRow(t_star,r_star,t,r);
                DrowMax = overlapMaxRow(t_star,r_star,t,r);
                DcolMin = overlapMinCol(t_star,r_star,t,r);
                DcolMax = overlapMaxCol(t_star,r_star,t,r);
                
                Drow = repmat(DrowMin:DrowMax,[1, DcolMax-DcolMin]);
                Dcol = repmat(DcolMin:DcolMax,[DrowMax-DrowMin, 1]);
                % Shift indices to center bi_star
                Rrow = mod(Drow + (row_star-1) - 1,m) + 1;
                Rcol = mod(Dcol(:) + (col_star-1) - 1,n) + 1; 
                
                % Zero overlapping atoms
                b(Rrow,Rcol,i,j) = 0; 
            end
        end
        fprintf('-] Spot placed\n')
        
        % Update max inner product
        [bi_star, i_star] = max(b(:));
        
        
        
    end
    
    % Update residual
    newR = x - forceMaskToZero(Ax_ft_2D(A0ft_stack,a),zMask);
    
    % Break if residual is not changing
    if (norm(R(:))-norm(newR(:)))/x_norm < params.delta
       break 
    end
    
    R = newR;
    
    fprintf('%3.3f\n', norm(R(:))/x_norm)
    k = k + 1;
end


end

