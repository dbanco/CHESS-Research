function [x_hat,nIter, timeSteps, errorSteps] = SolveFISTA_Circulant(A0ft_stack,b, varargin)
% b  - m x n polar ring image (required input)
% A0_stack - m x n x t x r unshifted gaussian basis matrices (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maxilambdam number of iterations
%         - DEFAULT 10000, if omitted or -1.
% lineSearchFlag - 1 if line search is to be done every iteration
%                - DEFAULT 0, if omitted or -1.
% continuationFlag - 1 if a continuation is to be done on the parameter lambda
%                  - DEFAULT 1, if omitted or -1.
% eta - line search parameter, should be in (0,1)
%     - ignored if lineSearchFlag is 0.
%     - DEFAULT 0.9, if omitted or -1.
% lambda - relaxation parameter
%    - ignored if continuationFlag is 1.
%    - DEFAULT 1e-3, if omitted or -1.
% outputFileName - Details of each iteration are dumped here, if provided.
%
% x_hat - estimate of coeeficient vector
% numIter - number of iterations until convergence

t0 = tic ;
STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_SUBGRADIENT;

stoppingCriterion = STOPPING_DEFAULT;
maxIter = 10000 ;
isNonnegative = 1;
tolerance = 1e-3;
[m,n,t,r] = size(A0ft_stack) ;
%x0 = ones(m,n,t,r) ;
xG = [];

%% Initializing optimization variables
t_k = 1 ; 
t_km1 = 1 ;
L0 = 1e6 ;
%G = A'*A ; Replace with two convolutions
nIter = 0 ;
c = AtR_ft_2D(A0ft_stack,b);
xk = ones(m,n,t,r) ;
% xk = A\b ; Might be good to init x this way if I can solve this quickly
L = L0 ;

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' lambdast be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'tolerance'
            tolerance = parameterValue;
        case 'lambda'
            lambda = parameterValue;
        case 'beta'
            beta = parameterValue;
        case 'maxiteration'
            maxIter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'initialization'
            xk = parameterValue;
            if ~all(size(xk)==[m,n,t,r])
                error('The dimension of the initial xk does not match.');
            end
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = nan(1,maxIter) ;
errorSteps = nan(1,maxIter) ;

if stoppingCriterion==STOPPING_GROUND_TRUTH && isempty(xG)
    error('The stopping criterion must provide the ground truth value of x.');
end

keep_going = 1 ;
nz_x = (abs(xk)> eps*10);
f = 0.5*norm(b-Ax_ft_2D(A0ft_stack,xk))^2 + lambda * norm(xk(:),1);
xkm1 = xk;
while keep_going && (nIter < maxIter)
    nIter = nIter + 1 ;
    
    % FISTA solution update
    yk = xk + ((t_km1-1)/t_k)*(xk-xkm1) ;
    
    % Enforce positivity
    if isNonnegative
        yk(yk<0) = 0;
    end
    
    % Compute gradient of f
    grad = AtR_ft_2D(A0ft_stack,Ax_ft_2D(A0ft_stack,yk)) - c ; % gradient of f at yk
    
    % Backtracking
    stop_backtrack = 0 ;
    while ~stop_backtrack
        
        gk = yk - (1/L)*grad ;
        if isNonnegative
            gk(gk<0) = 0;
        end
        xkp1 = soft(gk,lambda/L) ;
        
        % Compute objective at xkp1
        fit = Ax_ft_2D(A0ft_stack,xkp1);
        
        % Compute quadratic approximation at yk
        fit2 = Ax_ft_2D(A0ft_stack,yk);
        temp1 = norm(b(:)-fit(:))^2;
        temp2 = norm(b(:)-fit2(:))^2 +...
            (xkp1(:)-yk(:))'*grad(:) + (L/2)*norm(xkp1(:)-yk(:))^2 ;
        
        % Stop backtrack if objective <= quadratic approximation
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*beta ;
        end
        
    end

    timeSteps(nIter) = toc(t0) ;
    errorSteps(nIter) = norm(b(:)-fit(:))/norm(b(:)) ;
    
    disp(['Iter ', num2str(nIter),...
          '||x||_0 ', num2str(sum(abs(xkp1(:)) > 0)),...
          ' err ', num2str(errorSteps(nIter)) ]) ;
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            keep_going = norm(xG(:)-xkp1(:))>tolerance;
        case STOPPING_SUBGRADIENT
            sk = L*(yk-xkp1) +...
                 AtR_ft_2D(A0ft_stack,Ax_ft_2D(A0ft_stack,(xkp1-yk)));
            keep_going = norm(sk(:)) > tolerance*L*max(1,norm(xkp1(:)));
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            nz_x_prev = nz_x;
            nz_x = (abs(xkp1)>eps*10);
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                keep_going = (criterionActiveSet > tolerance);
            end
        case STOPPING_OBJECTIVE_VALUE
            % compute the stopping criterion based on the relative
            % variation of the objective function.
            prev_f = f;
            f = 0.5*norm(b-fit)^2 + lambda * norm(xk(:),1);
            criterionObjective = abs(f-prev_f)/(prev_f);
            keep_going =  (criterionObjective > tolerance);
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for PGBP.');
        case STOPPING_TIME
            keep_going = timeSteps(nIter) < maxTime ;
        otherwise
            error('Undefined stopping criterion.');
    end
    
    % lambda = max(eta*lambda,lambda_bar) ;
    
    % FISTA step size update 
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    t_km1 = t_k ;
    t_k = t_kp1 ;
    xkm1 = xk ;
    xk = xkp1 ;
end

x_hat = xk ;
timeSteps = timeSteps(1:nIter) ;
errorSteps = errorSteps(1:nIter) ;

function y = soft(x,T)
if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = sign(x).*y;
end