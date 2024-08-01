topDir = 'C:\Users\dpqb1\Documents\Outputs\';
outDir = 'mmpad_ring1_optFlow_OF_2_1_2024';
folderPath = fullfile(topDir,outDir);


files = dir(fullfile(folderPath, '*.mat'));


% Extract the file names and store them in a cell array
matFileNames = {files.name};

NN = numel(matFileNames);
rel_error = zeros(NN,1);
l1_norm = zeros(NN,1);
of_term = zeros(NN,1);
hs_term = zeros(NN,1);
lambdas_S = zeros(NN,1);
lambdas_OF = zeros(NN,1);
lambdas_HS = zeros(NN,1);
for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(folderPath,matFileNames{i}))

    D = outputs.D;
    N = outputs.N;
    M = outputs.M;
    T = outputs.T;
    y = outputs.y;
    X = outputs.X;
    K = outputs.K;
    opt = outputs.opt;
    scales = outputs.scales;
    lambda = outputs.lambda;
    lambda2 = outputs.lambda2;
    Uvel = outputs.Uvel;
    Vvel = outputs.Vvel;
    center = (M+1)/2;

    AD = reSampleCustomArrayCenter(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);
    err = sum((squeeze(y)-Yhat).^2,'all');
    % Compute error
%     rel_error(i) = sqrt(err/sum((squeeze(y)).^2,'all'))/T;
    rel_error(i) = err;
    % Comptue L1-norm
    l1_norm(i) = norm(X(:),1);
    [~,~,Fx,Fy,Ft] = computeHornSchunkDictPaperLS(X,K,Uvel,Vvel,opt.Smoothness/lambda2,opt.HSiters);
    [Jof, Jhs] = HSobjectivePaper(Fx,Fy,Ft,Uvel,Vvel,K,opt.Smoothness/lambda2);
    of_term(i) = Jof;
    hs_term(i) = Jhs;
    lambdas_S(i) = lambda;
    lambdas_HS(i) = opt.Smoothness;
    lambdas_OF(i) = lambda2;
end

%% Plot an L curve
close all
% Plot L-curve and selected parameter
        % [err_sort,ind] = sort(rel_error);
% l1_sort = l1_norm(ind);
% plot(l1_sort,err_sort,'o-')
% ylabel('||ADX - B||_2')
% xlabel('||X||_1')

x_data = of_term;
y_data = rel_error;
plot(x_data,y_data,'o-')
ylabel('Error')
xlabel('OF')

%% Select parameter
% criterion = (abs((err_sort-min(err_sort))/max(err_sort-min(err_sort)))) +...
%              0.5*(abs((l1_sort-min(l1_sort))/max(l1_sort-min(l1_sort))));

criterion = x_data/200 + y_data;

[minCrit, selInd] = min(criterion);

hold on
plot(x_data(selInd),y_data(selInd),'sr','MarkerSize',10)

% Get Lambda for chosen parameter
fprintf('lam_S = %0.2s \n',lambdas_S(selInd))
fprintf('lam_HS = %0.2s \n',lambdas_HS(selInd))
fprintf('lam_OF = %0.2s \n',lambdas_OF(selInd))

