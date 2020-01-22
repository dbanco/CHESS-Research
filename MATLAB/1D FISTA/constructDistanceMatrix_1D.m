function [ D ] = constructDistanceMatrix_1D( P, Threshold )
%constructDistanceMatrix_1D Summary of this function goes here
%   Detailed explanation goes here

% Construct distance matrix
N = P.num_var_t;
switch P.cost
    case 'l1'
        D = ones(N,N).*Threshold;
        for i = 1:P.num_var_t
            for ii=max([1 i-Threshold+1]):min([P.num_var_t i+Threshold-1])
                D(i,ii)= abs(i-ii);
            end
        end

        D = D./max(D(:));

    case 'l2'
        D = ones(N,N).*Threshold;
        for i = 1:P.num_var_t
            for ii=max([1 i-Threshold+1]):min([P.num_var_t i+Threshold-1])
                D(i,ii)= (i-ii)^2;
            end
        end

        D = D./max(D(:));

    case 'wass'
        D = ones(N,N).*Threshold;
        for i = 1:P.num_var_t
            for ii=max([1 i-Threshold+1]):min([P.num_var_t i+Threshold-1])
                D(i,ii)= P.var_theta(i) + P.var_theta(ii) -...
                    2*sqrt(P.var_theta(i)*P.var_theta(ii));
            end
        end
        D = D./max(D(:));

    case 'sqrt'
        D = ones(N,N).*Threshold;
        for i = 1:P.num_var_t
            for ii=max([1 i-Threshold+1]):min([P.num_var_t i+Threshold-1])
                D(i,ii)= sqrt(abs(i-ii));
            end
        end
        D = D./max(D(:));
end
end