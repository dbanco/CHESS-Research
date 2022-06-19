function [U, D1, D2, X1, X2, Vals] = cdc_ncp(Y,D1,D2,param)
%cdc_ncp Convolutional dictionary constrained nonnegative CP
     
% Set up status display for verbose operation
[N1,N2,N3] = size(Y);
Udict = {zeros(N1,param.R),zeros(N2,param.R),zeros(N3,param.R)};

count = 0;
Vals = zeros(1,11);
hstr = ['Itn   Obj    ncpErr    D1Err      D2Err     L3Err        '...
        'D1L1      D2L1  Iters:NCP, DL1, DL2  '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e    %4d %4d %4d \n';
nsep = 84;
if param.verbose && param.maxIters > 0
  disp(hstr);
  disp(char('-' * ones(1,nsep)));
end

switch param.Uinit
    case 'rand'
        U = cell(3,1);
        for i=1:3
            U{i}=rand(size(Y,i),param.R);
        end
    case 'ones'
        U = cell(3,1);
        for i=1:3
            U{i}=ones(size(Y,i),param.R);
    end
end

k = 1; noStop = 1;    
dlIter1 = param.opt1.MaxMainIter;
dlIter2 = param.opt2.MaxMainIter;
while (k <= param.maxIters) && noStop
    % Nonnegative CP
    if k == 1
        kappa = 0;
        phi1 = 0;
        phi2 = 0;
        ncpIters = 50;
    else
        kappa = param.kappa;
        phi1 = param.phi1;
        phi2 = param.phi2;
        ncpIters = param.ncpIters;
    end
    [U,ncpVals] = ConNCP(Y,param.R,...
                                'max_iter',ncpIters,...
                                'tol',param.ncpTol,...
                                'verbose',param.ncpVerbose,...
                                'Udict',Udict,...
                                'kappa',kappa,...
                                'phi1',phi1,...
                                'phi2',phi2,...
                                'init',U);           
    if k == 1
       param.opt1.MaxMainIter = 50; 
       param.opt2.MaxMainIter = 50; 
    else
       param.opt1.MaxMainIter = dlIter1; 
       param.opt2.MaxMainIter = dlIter2; 
    end
    % CBPDNDL Dictionary learning 
    [D1, X1, optinf1] = cbpdndl(D1,reshape(U{1},[1,size(U{1})]),param.lambda,param.opt1);
    Udict{1} = param.phi1*recon(D1,X1);
    [dl1Iters,D1Err, D1L1] = DLObj(optinf1,X1,Udict{1},U{1},param.phi1,param.lambda);

    [D2, X2, optinf2] = cbpdndl(D2,reshape(U{2},[1,size(U{2})]),param.lambda,param.opt2);
    Udict{2} = param.phi2*recon(D2,X2);
    [dl2Iters,D2Err, D2L1] = DLObj(optinf2,X2,Udict{2},U{2},param.phi2,param.lambda);

    F=ktensor(U);
    % Compute objective terms
    L3Err = param.kappa*0.5*norm(L(U{3}))^2;
    ncpErr = 0.5*norm(Y-full(F))^2;
    Obj = ncpErr + D1Err + D2Err + L3Err + D1L1 + D2L1;
    Vals(k,:) = [k, Obj, ncpErr, D1Err, D2Err, L3Err, D1L1, D2L1,...
                 ncpVals(end,1),dl1Iters,dl2Iters];
    
    % Update stopping criterion
    if k>1
        criterion = abs(Vals(k-1,2)-Vals(k,2))/Vals(k-1,2);
    else
        criterion = 1;
    end     
    % Check for convergence
    if  ((criterion >= 0) &&...
        (criterion < param.stopCrit)) || isnan(criterion)
        count = count + 1;
        if count >= 5
            noStop = 0;
        end
    else
        count = 0;
    end
    if param.verbose
        fprintf(sfms, Vals(k,:));
    end
    % Carry over ADMM vars
    param.opt1 = carryOverOpt(param.opt1,optinf1);
    param.opt2 = carryOverOpt(param.opt2,optinf2);

    k = k + 1;
end
if param.verbose
    disp(['----end--',char('-' * ones(1,nsep))]);
end

function [Iters,Err,L1] = DLObj(optinf,X,Uhat,U,phi,lambda)
Iters = size(optinf.itstat,1);
Err = 0.5*phi*norm(Uhat-U)^2;
L1 = phi*lambda*abs(sum(X(:)));
return

function opt = carryOverOpt(opt,optinf)
opt.Y0 = optinf.Y; 
opt.sigma = optinf.sigma;
return

function Y = recon(D,X)
Df = fft2(D,1,size(X,2));
Xf = fft2(X);
Yf = sum(bsxfun(@times,Df,Xf),3);
Y = squeeze(real(ifft2(Yf)));
return

function Y = L(X)
[T,R] = size(X);
Y = zeros(T,R);
Y(1,:) = -2*X(1,:) + 2*X(2,:);
Y(T,:) = 2*X(T-1,:) - 2*X(T,:);
for t = 2:(T-1)
    Y(t,:) = -2*X(t,:) + X(t+1,:) + X(t-1,:);
end
return