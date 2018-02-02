function [x_hat_array] = cyclic_optimization(A0ft_stack,b_array,params)
%cyclic_optimization Alternately uses FISTA algorithm across HEXD data array 
%
% Inputs:
% b_array    - (N x M) Cell array of(m x n) polar ring images
% A0ft_stack - (m x n x t x r) fft2 of unshifted gaussian basis matrices
% params     - struct containing the following field
%   lambda - l1 penalty parameter > 0
%   gamma - spatial smoothness regularization parameter > 0
%   L - initial Lipschitz constant > 0
%   beta - backtracking parameter > 1
%   stoppingCriterion - integer indicated stopping criterion (1,2,3)
%   tolerance - tolerance for stopping criterion
%   maxIter - maximum number of iterations
%   maxCycles - maximum number of out iterations
%   isNonnegative - flag to enforce nonnegative solution
%   x_init - initial guess of solution
%
% Outputs:
% x_hat_array - (m x n x t x r) solution 
[m,n,t,r] = size(A0ft_stack);
[N,M] = size(b_array);
x_hat_array = cell(N,M);
error_array = zeros(N,M);

% Initial cycle without regularization (All in parallel)
    for i = 1:N
        for j = 1:M
            [x_hat,err,~,l0] = FISTA_Circulant(A0ft_stack,...
                                            b_array{i,j},...
                                            params);
            x_hat_array{i,j} = x_hat;
            error_array(i,j) = err(end)
            l0_array = l0(end)
        end
    end

% Cycles with regularization
for k = 1:params.maxCycles
    for i = 1:N
        for j = 1:M
            fprintf('Cycle: %i;  Point %i,%i',k,i,j)
            x_neighbors = neighbors_vdf(x_hat_array,i,j);
            params.x_init = x_hat_array{i,j};
            [x_hat,err,~,l0] = space_FISTA_Circulant(A0ft_stack,...
                                                  b_array{i,j},...
                                                  x_neighbors,params);
            x_hat_array{i,j} = x_hat;
            error_array(i,j) = err(end)
            l0_array = l0(end)
        end
    end
end
