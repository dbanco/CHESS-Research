function a = GCMP(A0ft_stack, x, params)
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
    
    figure(2)
    [ha2,pos2] = tight_subplot(4,1,[0.01,0.01],[0.01,0.01],[0.01,0.01])
    axes(ha2(1))
    imshow(x,'DisplayRange',[0 uplim],'Colormap',jet)
    p_ind = 2;
    pause(0.05)
end

% Precompute norm of data
x_norm = norm(x(:));

num_spots = 0;
while (norm(R(:))/x_norm) > params.epsilon
    % Compute inner products
    b = AtR_ft_2D(A0ft_stack,R);
    
    % Enforce nonnegativity
    if params.isNonnegative
        b(b<1e-12) = 0;
    end
    
    % Find max inner product
    [bi_star, i_star] = max(b(:));
    
%     while bi_star > 1e-1

        % Update solution
        a(i_star) = a(i_star) + bi_star;
        
        num_spots = num_spots + 1;
        
        if params.showImage
            fit = forceMaskToZero(Ax_ft_2D(A0ft_stack,a),zMask);
            figure(1)
            subplot(2,1,2)
            imshow(fit,'DisplayRange',[0 uplim],'Colormap',jet)
            title('Da')
            pause(0.05)
        end
        
        if params.showImage && sum(num_spots==[8,20,45])
            fit = forceMaskToZero(Ax_ft_2D(A0ft_stack,a),zMask);
            figure(2)
            axes(ha2(p_ind))
            imshow(fit,'DisplayRange',[0 uplim],'Colormap',jet)
            title([nums2str(num_spots)])
            p_ind = p_ind + 1;
        end

        % Convert linear index i_star to 4D dictionary index
        [row_star,col_star,t_star,r_star] = ind2sub([m,n,t,r],i_star);
        
        % Loop over all dictionary atom shaps to remove overlapping atoms
        fprintf('[-')
%         for i = 1:t
%             fprintf('-%i-',round(100*i/t))
%             for j = 1:r
%                 % Convolve atoms to identify overlap
%                 atom_ft_star = squeeze(A0ft_stack(:,:,t_star,r_star));
%                 atom_ft_ij = squeeze(A0ft_stack(:,:,t,r));
%                 conv_atom = real(ifft2(atom_ft_star.*atom_ft_ij));
%                 [Drow,Dcol] = find(conv_atom > 2e-1); 
%            
%                 
%                 % Shift indices to center bi_star
%                 Rrow = mod(Drow + (row_star),m) + 1;
%                 Rcol = mod(Dcol + (col_star),n) + 1; 
%                 
%                 Ind = sub2ind([m,n],Rrow,Rcol);
%                 aaa = zeros(m,n);
%                 aaa(Ind) = 1;
%                 aaa(row_star,col_star) = 0;
%                 
%                 % Zero overlapping atoms
%                 b_slice = b(:,:,i,j);
%                 b_slice(Ind) = 0;
%                 b(:,:,i,j) = b_slice; 
%                 
%                 figure(11)
%                 subplot(3,1,1)
%                 imagesc(aaa)
%                 subplot(3,1,2)
%                 imagesc(conv_atom > 2e-1)
%                 subplot(3,1,3)
%                 imagesc(b(:,:,i,j))
%             end
%         end
%         fprintf('-] Spot placed\n')
%         
%         % Update max inner product
%         [bi_star, i_star] = max(b(:));
%         
%     end
    
    % Update residual
    newR = x - forceMaskToZero(Ax_ft_2D(A0ft_stack,a),zMask);
    
    % Break if residual is not changing
    if (norm(R(:))-norm(newR(:)))/x_norm < params.delta
       break 
    end
    
    R = newR;
    
    fprintf('%3.3f\n', (norm(R(:))/x_norm))
    k = k + 1;
end


end

