function x = fista_circulant_2D(A0ft_stack, b, x, L,l1_ratio, maxit, eps, positive, verbose, benchmark)
%fista_circulant_2D Summary of this function goes here
%   FISTA algorithm using circulant matrix-vector product subroutines
    % A0 is a bunch of slices indexed by variance and radius
if nargin < 9
    benchmark = 0;
end
if nargin < 8
    verbose = 0;
end
if nargin < 7
    positive = 1;
end
if nargin < 6
    eps = 10^(-7);
end
if nargin < 5
    maxit = 800;
end

if benchmark 
    tic;
end

Linv = 1/L;

if x == 0
    x = zeros(size(A0ft_stack));
end
t = 1;
z = x;

for it = 1:maxit
    xold = x;

    % Arrange x coefficents as matrix in fourier domain 
    R = b - Ax_ft_2D(A0ft_stack,z);
    z = z + AtR_ft_2D(A0ft_stack,R)*Linv;

    % Enforce positivity on coefficients
    if positive
        z(z<0 ) = 0;
    end
    x = soft_thresh(z,l1_ratio*Linv);

    t0 = t;
    t = (1 + sqrt(1 + 4*t^2))/ 2;
    z = x + ((t0 - 1)/t)*(x - xold);

    if positive
        z(z<0) = 0;
    end

    criterion = sum(abs(x(:) - xold(:)))/numel(x);

    if verbose == 1
        disp(['Iteration  ' + num2str(it) +' of '+ num2str(maxit)])
    end
    if verbose == 2
        if (mod(it, 10) == 1)
            fit = Ax_ft_2D(A0ft_stack,x);
            sqrt_res = norm(b(:) - fit(:));
            L0 = sum(x(:)>0);
            disp(['Iter ',num2str(it), ' of ', num2str(maxit), ...
                  ', Err ', num2str(sqrt_res/norm(b(:)))...
                  ', ||x||_0 ', num2str(L1)])
        end
    end
    if verbose == 3
        fit = Ax_ft_2D(A0ft_stack,x);
        sqrt_res = norm(b(:) - fit(:));
        L0 = sum(x(:)>0);
        disp(['Iter ',num2str(it), ' of ', num2str(maxit), ...
              ', Err ', num2str(sqrt_res/norm(b(:)))...
              ', ||x||_0 ', num2str(L0)])
    end
    if(criterion < eps)
        break
    end

    if benchmark 
        total_time = toc;
    end

end

function y = soft_thresh(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end