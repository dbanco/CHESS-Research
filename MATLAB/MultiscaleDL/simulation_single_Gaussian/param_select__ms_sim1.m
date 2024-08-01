topDir = 'C:\Users\dpqb1\Documents\Outputs2024\';
outDir = 'gaus_example_1_29_24_X0_D0_V00_sig_2\';
folderPath = fullfile(topDir,outDir);


files = dir(fullfile(folderPath, '*.mat'));


% Extract the file names and store them in a cell array
matFileNames = {files.name};

NN = numel(matFileNames);
rel_error = zeros(NN,1);
l1_norm = zeros(NN,1);
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
%     rel_error(i) = sqrt(err/sum((squeeze(y)).^2,'all'))/T;
    rel_error(i) = sqrt(err);
    % Comptue L1-norm
    l1_norm(i) = norm(X(:),1);
    
end

% Plot L-curve and selected parameter
[err_sort,ind] = sort(rel_error);
l1_sort = l1_norm(ind);
loglog(l1_sort(1:end-3),err_sort(1:end-3),'o-')
ylabel('Error')
xlabel('l_1-norm')

criterion = abs((err_sort-min(err_sort))/max(err_sort-min(err_sort))) +...
            abs((l1_sort-min(l1_sort))/max(l1_sort-min(l1_sort)));
[minCrit, selInd] = min(criterion);

hold on
loglog(l1_sort(3),err_sort(3),'sr','MarkerSize',10)

% Get Lambda for chosen parameter
load(fullfile(folderPath,matFileNames{ind(selInd)}))
outputs.lambda

