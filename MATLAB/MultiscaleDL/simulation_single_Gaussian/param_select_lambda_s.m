function [lambda_s,selInd] = param_select_lambda_s(outputDir,tradeoff,scaleP)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

files = dir(fullfile(outputDir, '*.mat'));

% Extract the file names and store them in a cell array
matFileNames = {files.name};

NN = numel(matFileNames);
rel_error = zeros(NN,1);
l1_norm = zeros(NN,1);
lambda_s_vec = zeros(NN,1);
for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(folderPath,matFileNames{i}))

    D = outputs.D;
    N = outputs.N;
    M = outputs.M;
    T = outputs.T;
    y = outputs.y;
    X = outputs.X;
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
    % Get lambda parameter
    lambda_s_vec(i) = outputs.lambda;
end

[lambda_s_sort,ind] = sort(lambda_s_vec);
err_sort = rel_error(ind);
l1_sort = l1_norm(ind);

% Normalized origin distance criterion
criterion = tradeoff*abs((err_sort-scaleP(1))/scaleP(2)) +...
                abs((l1_sort-scaleP(3))/scaleP(4));
[~, selInd] = min(criterion);
lambda_s = lambda_s_sort(selInd);

end