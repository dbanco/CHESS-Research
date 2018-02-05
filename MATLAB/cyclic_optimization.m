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

[N,M] = size(b_array);
x_hat_array = cell(N,M);
error_array = zeros(N,M);

% Set up indices
all = 1:N*M;
odd = 1:2:N*M;
even = 2:2:N*M;
odd_x_hat_array = cell(numel(odd),1);
even_x_hat_array = cell(numel(even),1);
odd_error_array = zeros(numel(odd),1);
even_error_array = zeros(numel(even),1);
odd_b_array = b_array(odd);
even_b_array = b_array(even);
% Initial cycle without regularization (All in parallel)
parfor idx = all
    [x_hat,err,~,l0] = FISTA_Circulant(A0ft_stack,...
                                    b_array{idx},...
                                    params);
    x_hat_array(idx) = {x_hat};
    error_array(idx) = err(end);
    l0_array = l0(end)
end

% Cycles with regularization
for k = 1:params.maxCycles
    parfor idx = 1:numel(odd)
        [i,j] = ind2sub([N,M],odd(idx));
        fprintf('Cycle: %i;  Point %i,%i',k,i,j)
        x_neighbors = neighbors_vdf(x_hat_array,i,j);
        local_params = params;
        local_params.x_init = x_hat_array{i,j};
        [x_hat,err,~,l0] = space_FISTA_Circulant(A0ft_stack,...
                                              odd_b_array{idx},...
                                              x_neighbors,local_params);
        odd_x_hat_array(idx) = {x_hat};
        odd_error_array(idx) = err(end);
        l0(end)
   
    end
    x_hat_array = fill_array(x_hat_array,odd_x_hat_array,odd);
    error_array = fill_array(error_array,odd_error_array,odd);
    parfor idx = 1:numel(even)
        [i,j] = ind2sub([N,M],even(idx));
        fprintf('Cycle: %i;  Point %i,%i',k,i,j)
        x_neighbors = neighbors_vdf(x_hat_array,i,j);
        local_params = params;
        local_params.x_init = x_hat_array{i,j};
        [x_hat,err,~,l0] = space_FISTA_Circulant(A0ft_stack,...
                                              even_b_array{idx},...
                                              x_neighbors,params);
        even_x_hat_array(idx) = {x_hat};
        even_error_array(idx) = err(end);
        l0(end)
    end
    x_hat_array = fill_array(x_hat_array,even_x_hat_array,even);
    error_array = fill_array(error_array,even_error_array,even);
end
