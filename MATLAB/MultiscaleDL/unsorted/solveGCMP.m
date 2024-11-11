function x = solveGCMP(D, b, opt)
%GCMP Group Convolutional Matching Pursuit
% Solves a relaxed version of the convolutional sparse coding problem 
%    min ||x||_0,inf
%    s.t. b = Dx 
% referred to as the P^epsilon_0,inf problem
%    min ||x||_0,inf
%    s.t. ||b - Dx||^2 <= epsilon
% where generally speaking,
% D is a BCCB dictionary (block circulant with circulant blocks)
% a are linear coefficients defing the solution
% b is input data to be coded
% epsilon is an error bound.
%
% Implemented mostly as described in
% A Greedy Approach to l_0, Based Convolutional Sparse Coding
% by Elad Plaut and Raja Giryes
%
% Inputs:
% b - (N1 x N2 X T) data
% D - (M1 x M2 x K) Dictionary
% params - struct containing the following field
%   epsilon - relative error stopping parameter in (0,1)
%   delta - change in relative error stopping parameter in (0,1)
%   zeroMask - (#pixels x 2) indices of unobserved pixels 
%   isNonnegative - flag to enforce nonnegative solution
%   showImage - flag to display image as each spot is placed
%
% Outputs:
% x - (N1 x N2 X K X T) solution 
[~,~,KJ] = size(D);
[N1,N2,~,T] = size(b);
bpad = padarray(b,[0 N2-1 0 0],0,'pre');
[~,N2pad,~] = size(bpad);

x = zeros(N1,N2pad,KJ,T);

wrapN = @(i, N) (1 + mod(i-1, N));

dictFilterSizes = zeros(2,KJ);
dictFilterSizes(1,:) = N1;

scales = genRationals([0;1],[1;1],16,16, 1/8);
for i = 1:KJ
    dictFilterSizes(2,i) = floor(N2*scales(1,i)/scales(2,i));
end

Df = fft2(D);
bf = fft2(bpad);
Rf = bf;

layers = 0;
num_spots = 0;
while (norm(Rf(:)) > opt.epsilon) && (layers < opt.maxLayers)
    % Compute inner products
    y = ifft2(bsxfun(@times,conj(Df),Rf),'symmetric');
    y(:,1:N2-1,:,:) = 0;
    
    % Enforce nonnegativity
    if opt.NonNegCoef
        y(y<1e-12) = 0;
    end
    
    % Find max inner product
    [yi_star, i_star] = max(y,[],[1,2,3],'linear');

    while max(y(:)) > opt.epsilon  

        % Update solution
        x(i_star) = x(i_star) + yi_star;
            
        num_spots = num_spots + 1;
        
        for t = 1:T
            % Convert linear indices i_star to 4D dictionary index
            [n1_star,n2_star,k_star,t_star] = ind2sub([N1,N2pad,KJ,T],squeeze(i_star(t)));
    

            % Remove spots overlapping
            M1star = dictFilterSizes(1,k_star);
            M2star = dictFilterSizes(2,k_star);
            for k = 1:KJ
                Mk = dictFilterSizes(2,k);
                low1 = n1_star - Mk + 1;
                low1 = max(1,low1);
                hi1 = n1_star + M1star' -1;
                hi1 = min(hi1,N1);
                low2 = n2_star - Mk + 1;
                low2 = max(1,low2);
                hi2 = n2_star + M2star' -1;
                hi2 = min(N2pad,hi2);
                
                y(low1:hi1,low2:hi2,k,t) = 0;
            end
        end
                    
        % Update max inner product
        [yi_star, i_star] = max(y,[],[1,2,3],'linear');

    end
    
    % Update residual
    bhatf = sum(bsxfun(@times,Df,fft2(x)),3);
    bhatf(:,1:N2-1,:,:) = 0;
    newRf = bf - bhatf;
    
    % Break if residual is not changing
    if norm(Rf(:))-norm(newRf(:)) < opt.delta
       break 
    end

    Rf = newRf;
    
    fprintf('%3.3f\n', norm(Rf(:))  )
    layers = layers + 1;

end


end

