function [lambda_of,outInd] = param_select_lambda_of(outputDir,tradeoff,scaleP,fig_num)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

files = dir(fullfile(outputDir, '*.mat'));

% Extract the file names and store them in a cell array
matFileNames = {files.name};

NN = numel(matFileNames);
rel_error = zeros(NN,1);
l1_norm = zeros(NN,1);
of_obj = zeros(NN,1);
hs_obj = zeros(NN,1);
lambda_of_vec = zeros(NN,1);
for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(outputDir,matFileNames{i}))

    D = outputs.D;
    N = outputs.N;
    M = outputs.M;
    T = outputs.T;
    y = outputs.y;
    X = outputs.X;
    K = outputs.K;
    scales = outputs.scales;
    center = (M+1)/2;

    AD = reSampleCustomArrayCenter(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);
    err = sum((squeeze(y)-Yhat).^2,'all');
    % Compute error
    rel_error(i) = sqrt(err);
    % Compute L1-norm
    l1_norm(i) = norm(X(:),1);
    % Compute OF objective
    Xpad = padarray(squeeze(outputs.X),[1 1 1],0,'pre');
    Fx = diffxHS(Xpad);
    Fy = diffyHS(Xpad);
    Ft = difftHS(Xpad);
    [Jof, Jhs] = HSobjectivePaper(Fx,Fy,Ft,outputs.Uvel,outputs.Vvel,K,outputs.opt.Smoothness/outputs.lambda2);
    of_obj(i) = Jof;
    hs_obj(i) = Jhs;
    % Get lambda parameter
    lambda_of_vec(i) = outputs.lambda2;
end

[lambda_of_sort,ind] = sort(lambda_of_vec);
err_sort = rel_error(ind);
l1_sort = l1_norm(ind);
of_sort = of_obj(ind);

% Normalized origin distance criterion
criterion = tradeoff*(abs((err_sort-scaleP(1))/scaleP(2))+abs((l1_sort-scaleP(3))/scaleP(4))) +...
                abs((of_sort-scaleP(5))/scaleP(6));
[~, selInd] = min(criterion);
lambda_of = lambda_of_sort(selInd);

if nargin > 3
    figure(fig_num)
    loglog(of_sort(1:end),err_sort(1:end),'o-')
    ylabel('Error')
    xlabel('OF Obj')
    hold on
    loglog(of_sort(selInd),err_sort(selInd),'sr','MarkerSize',10)
end

outInd = ind(selInd);

end