function x = solveGCMP(Df, y, opt)
%GCMP Group Convolutional Matching Pursuit
% Solves a relaxed version of the convolutional sparse coding problem 
%    min ||x||_0,inf
%    s.t. y = Dx 
% referred to as the P^epsilon_0,inf problem
%    min ||x||_0,inf
%    s.t. ||y - Dx||^2 <= epsilon
% where generally speaking,
% Df is contains the fft2 of dictionary atoms
% x are linear coefficients defing the solution
% y is input data to be coded
%
% Implemented mostly as described in
% A Greedy Approach to l_0, Based Convolutional Sparse Coding
% by Elad Plaut and Raja Giryes
%
% Inputs:
% y - (N1 x N2) data
% Df - (N1 x N2 x K) fft of Dictionary
% opt - struct containing the following fields
%   maxLayers - maximum number of nonoverlapping layers to place (integer)
%   epsilon - relative error stopping parameter in (0,1)
%   delta - change in relative error stopping parameter in (0,1)
%   NonNegCoef - flag to enforce nonnegative solution (0 or 1)
%   DictFilterSizes - (2 x K) array containing the widths of dictionary entries so
%   that we can check for overlap 
% Outputs:
% x - (N1 x N2 X K) solution 

[~,~,K] = size(Df);
[N1,N2] = size(y);
x = zeros(N1,N2,K);

wrapN = @(i, N) (1 + mod(i-1, N));

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
    
    % Continue until no more large enough nonoverlapping functions can be placed 
    while max(b(:)) > 0.25*norm(y(:)) 

        % Update solution
        x(i_star) = x(i_star) + bi_star;
            
        num_spots = num_spots + 1;
        
        % Convert linear indices i_star to 4D dictionary index
        [n1_star,n2_star,k_star,~] = ind2sub([N1,N2,K],squeeze(i_star));
        hold on
        plot(n2_star,n1_star,'o')

        % Remove spots overlapping
        M1star = opt.DictFilterSizes(1,k_star);
        M2star = opt.DictFilterSizes(2,k_star);
        sum(abs(b(:))>0)
        for k = 1:K
            Mk1 = opt.DictFilterSizes(1,k) + M1star;
            Mk2 = opt.DictFilterSizes(2,k) + M2star;
            low1 = n1_star - Mk1 + 1;
            hi1 = n1_star + Mk1 -1;
            low2 = n2_star - Mk2 + 1;
            hi2 = n2_star + Mk2 -1;
            b(wrapN(low1:hi1,N1),wrapN(low2:hi2,N2),k) = 0;

        end
        sum(abs(b(:))>0)
                
        % Update max inner product
        [bi_star, i_star] = max(b,[],[1,2,3],'linear');
    
    end
    
    % Update residual
    yfhat = sum(bsxfun(@times,Df,fft2(x)),3);
    newRf = yf - yfhat;
    
    % Break if residual is not changing
    if norm(Rf(:))-norm(newRf(:)) < opt.delta
       break 
    end

    Rf = newRf;
    
    fprintf('%3.3f\n', norm(Rf(:))  )
    layers = layers + 1;

end


end

