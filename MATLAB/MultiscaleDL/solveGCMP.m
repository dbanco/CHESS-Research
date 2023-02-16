function x = solveGCMP(D, y, opt)
%GCMP Group Convolutional Matching Pursuit
% Solves a relaxed version of the convolutional sparse coding problem 
%    min ||x||_0,inf
%    s.t. y = Dx 
% referred to as the P^epsilon_0,inf problem
%    min ||x||_0,inf
%    s.t. ||y - Dx||^2 <= epsilon
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
% Inputs:
% y - (N1 x N2 X T) data
% D - (M1 x M2 x K) Dictionary
% params - struct containing the following field
%   epsilon - relative error stopping parameter in (0,1)
%   delta - change in relative error stopping parameter in (0,1)
%   zeroMask - (#pixels x 2) indices of unobserved pixels 
%   isNonnegative - flag to enforce nonnegative solution
%   showImage - flag to display image as each spot is placed
%
% Outputs:
% x - (N1 x N2 X T X K) solution 
[~,~,K] = size(D);
[N1,N2,T] = size(y);
x = zeros(N1,N2,K,T);

wrapN = @(i, N) (1 + mod(i-1, N));

Df = fft2(D);
yf = fft2(y);
Rf = yf;

layers = 0;
num_spots = 0;
while (norm(Rf(:)) > opt.epsilon) && (layers < opt.maxLayers)
    % Compute inner products

    b = ifft2(bsxfun(@times,conj(Df),Rf),'symmetric');
    
    % Enforce nonnegativity
    if opt.NonNegCoef
        b(b<1e-12) = 0;
    end

    % Find max inner product
    [bi_star, i_star] = max(b,[],[1,2,3],'linear');
    
    while max(b(:)) > 1e-8        

        % Update solution
        x(i_star) = x(i_star) + bi_star;
            
        num_spots = num_spots + 1;
        
        % Convert linear indices i_star to 4D dictionary index
        [n1_star,n2_star,k_star,~] = ind2sub([N1,N2,K,T],squeeze(i_star));
    

        % Remove spots overlapping
        M1star = opt.DictFilterSizes(1,k_star);
        M2star = opt.DictFilterSizes(2,k_star);
        for k = 1:K
            Mk = opt.DictFilterSizes(2,k);
            low1 = n1_star - Mk + 1;
            hi1 = n1_star + M1star' -1;
            low2 = n2_star - Mk + 1;
            hi2 = n2_star + M2star' -1;
            for t = 1:T
                b(wrapN(low1(t):hi1(t),N1),wrapN(low2(t):hi2(t),N2),k,t) = 0;
            end
        end
                
        % Update max inner product
        [bi_star, i_star] = max(b,[],[1,2,3],'linear');

    end
    
    % Update residual
    newRf = yf - sum(bsxfun(@times,Df,fft2(x)),3);
    
    % Break if residual is not changing
    if norm(Rf(:))-norm(newRf(:)) < opt.delta
       break 
    end

    Rf = newRf;
    
    fprintf('%3.3f\n', norm(Rf(:))  )
    layers = layers + 1;

end


end

