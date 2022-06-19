function [U,Vals] = ConNCP(Y,R,varargin)
%ConNCP Summary of this function goes here
%   Detailed explanation goes here
% set parameters
params = inputParser;
params.addParameter('tol'              ,1e-4       ,@(x) isscalar(x) & x > 0 );
params.addParameter('max_iter'         ,200        ,@(x) isscalar(x) & x > 0 );
params.addParameter('init'             ,cell(0)    ,@(x) iscell(x) );
params.addParameter('verbose'          ,0          ,@(x) isscalar(x) & x >= 0 );
params.addParameter('orderWays'        ,[]);
params.addParameter('Udict'            ,cell(0)    ,@(x) iscell(x) );
params.addParameter('kappa'            ,1);
params.addParameter('phi1'            ,1);
params.addParameter('phi2'            ,1);

params.parse(varargin{:});

% copy from params object
P = params.Results;
P.nWay = ndims(Y);
P.R = R;
P.size = size(Y);

if isempty(P.orderWays)
    P.orderWays = [1:P.nWay]; 
end

Vals = zeros(1,6);
hstr = ['Itn   Obj    ncpErr    D1Err      D2Err     L3Err  '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e \n';
nsep = 84;
if P.verbose && P.max_iter > 0
  disp(hstr);
  disp(char('-' * ones(1,nsep)));
end

%Initialize
if ~isempty(P.init)
    U = P.init;
    P.init_type = 'User provided';
    P.init = cell(0);
else
    Finit = cell(P.nWay,1);
    for i=1:P.nWay
        Finit{i}=rand(size(Y,i),R);
    end
    U = Finit;
    P.init_type = 'Randomly generated';
end

d = ones(1,P.R);
for k=1:P.nWay-1
    curWay = P.orderWays(k);
    Fnorm2 = sqrt(sum(U{curWay}.^2,1));
    U{curWay} = U{curWay}./repmat(Fnorm2,size(U{curWay},1),1);
    d = d .* Fnorm2;
end
curWay = P.orderWays(end);
U{curWay} = U{curWay}.*repmat(d,size(U{curWay},1),1);

cellFF = cell(P.nWay,1);
for k=1:P.nWay
    cellFF{k} = U{k}'*U{k};
end

F_kten = ktensor(U);
Y_ten = tensor(Y);
% main iterations
cntu = 1;
epsilon = 1e-16;
for iter=1:P.max_iter
    % Init the scaling factors (gamma)
    d = sum(U{P.orderWays(end)}.^2,1);
    P.Udict{3} = P.kappa*SS(U{3});
    
    for k=1:P.nWay
        curWay = P.orderWays(k);
        ways = 1:P.nWay;
        ways(curWay)='';
        % Calculate Fnew = Y_(n) * khatrirao(all U except n, 'r').
        YF = mttkrp(Y_ten,U,curWay);
        % Compute the inner-product matrix
        FF = ones(P.R,P.R);
        for i = ways
            FF = FF .* cellFF{i};
        end
        if k < P.nWay
            for j = 1:P.R
                U{curWay}(:,j) = max( P.Udict{curWay}(:,j) + d(j)*U{curWay}(:,j) + YF(:,j) - U{curWay} * FF(:,j),epsilon);
                U{curWay}(:,j) = U{curWay}(:,j) ./ norm(U{curWay}(:,j)); 
            end
        else
            for j = 1:P.R
                U{curWay}(:,j) = max( P.Udict{curWay}(:,j) + U{curWay}(:,j) + YF(:,j) - U{curWay} * FF(:,j),epsilon)/(1+P.kappa);
%             U{curWay}(:,j) = max( P.Udict{curWay}(:,j) + U{curWay}(:,j) + XF(:,j) - U{curWay} * FF(:,j),epsilon);
            end
        end
        cellFF{curWay} = U{curWay}'*U{curWay};
    end
         
    % Compute objective terms
    D1Err =  P.phi1*0.5*norm(U{1}-P.Udict{1})^2;
    D2Err =  P.phi2*0.5*norm(U{2}-P.Udict{2})^2;
    L3Err = P.kappa*0.5*norm(U{3}-SS(F_kten.U{3}))^2;
    F_kten = ktensor(U);
    ncpErr = 0.5*max( norm(full(F_kten)-Y)^2 ,0);
    Obj = ncpErr + D1Err + D2Err + L3Err;
    Vals(iter,:) = [iter, Obj, ncpErr, D1Err, D2Err, L3Err];

    % Update stopping criterion
    if iter>1
        criterion = abs(Vals(iter-1,2)-Vals(iter,2))/Vals(iter-1,2);
    else
        criterion = 1;
    end
    % Check for convergence
    if  criterion < P.tol
        count = count + 1;
        if count >= 5
            cntu = 0;
        end
    else
        count = 0;
    end
    if P.verbose
        fprintf(sfms, Vals(iter,:));
    end

    if cntu==0, break; end
 
end

end

function Y = SS(X)
[T,R] = size(X);
Y = zeros(T,R);
Y(1,:) = X(2,:);
Y(T,:) = X(T-1,:);
for i = 2:T-1
    Y(i,:) = 0.5*(X(i-1,:) + X(i+1,:));
end
end