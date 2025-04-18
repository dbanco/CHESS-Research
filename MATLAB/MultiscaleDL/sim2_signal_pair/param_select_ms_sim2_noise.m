sigmas = 0:0.01:0.1;
NN = numel(sigmas);
selected_lam = zeros(NN,1);
selected_lam_true = zeros(NN,1);
close all
for n = 2:NN
    topDir = 'C:\Users\dpqb1\Documents\Outputs2024_8_29\';
    outDir = "signal_pair_8_29_24_X0_D0_V00_sig_"+num2str(n);
    folderPath = fullfile(topDir,outDir);
    
    files = dir(fullfile(folderPath, '*.mat'));
    [y,y_true,N,M,T] = generate_signal_pair(sigmas(n));
    % Extract the file names and store them in a cell array
    matFileNames = {files.name};
    
    NN = numel(matFileNames);
    error = zeros(NN,1);
    error_true = zeros(NN,1);
    l1_norm = zeros(NN,1);
    lambda_s_vec = zeros(NN,1);
    lambda_ind_vec = zeros(NN,1);
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
        err_true = sum((squeeze(y_true)-Yhat).^2,'all');
        % Compute error
        error(i) = sqrt(err);
        error_true(i) = sqrt(err_true);
        % Compute L1-norm
        l1_norm(i) = norm(X(:),1);
        % Get lambda parameter
        lambda_s_vec(i) = outputs.lambda;    
    end
    
    % Plot L-curve and selected parameter
    [lambda_s_sort,ind] = sort(lambda_s_vec);
    err_sort = error(ind);
    err_true_sort = error_true(ind);
    l1_sort = l1_norm(ind);
    fig_ps = figure(1);
    plot(l1_sort(1:end),err_sort(1:end),'o-')
    ylabel('Error')
    xlabel('l_1-norm')

    % Normalized origin distance criterion
    if n == 2
        minErr = min(err_sort);
        maxErr = max(err_sort-minErr);
        minL1 = min(l1_sort);
        maxL1 = max(l1_sort-minL1);
    end

    criterion = 1.5*abs((err_sort-minErr)/maxErr) +...
                abs((l1_sort-minL1)/maxL1);
    [minCrit, selInd] = min(criterion);
    
    [~,selInd2] = min(error_true);

    hold on
    plot(l1_sort(selInd),err_sort(selInd),'sr','MarkerSize',10);
    saveas(fig_ps,fullfile(folderPath,'param_select.png'));

    % Get Lambda for chosen parameter
    selected_lam(n) = lambda_s_sort(selInd);
    selected_lam_true(n) = lambda_s_sort(selInd2);
% 
%     figure
%     plot(curvature)

end
selected_lam
selected_lam_true

%% Next plot the true recovery as a function of SNR
sigmas = 0:0.01:0.1;
meanSNR = zeros(numel(sigmas),1);
for n = 2:numel(sigmas)
    [y,y_true,N,M,T] = generate_signal_pair(sigmas(n));
    
    SNR = zeros(T,1);
    nPwr = zeros(T,1);
    sigPwr = zeros(T,1);
    for t = 1:T
        SNR(t) = norm(y_true(:,t))/norm(y(:,t)-y_true(:,t));
        nPwr(t) = norm(y(:,t)-y_true(:,t));
        sigPwr(t) = norm(y_true(:,t));
    end
    
    meanSNR(n) = mean(SNR);
end
%%
NN = numel(sigmas);
trueErr = zeros(NN,1);
dataErr = zeros(NN,1);
noiseNorm = zeros(NN,1);
noiseNorm2 = zeros(NN,1);
outDir = "C:\Users\dpqb1\Documents\Outputs2024_8_29\";

for n = 2:NN
    outputDir = "signal_pair_8_29_24_X0_D0_V00_sig_"+num2str(n);
    dataFile = sprintf("output_j1_sig_%1.2s_lam1_%1.2s_lam2_0.00e+00",sigmas(n),selected_lam(n));
    [~,y_true,~,~,T] = generate_signal_pair(sigmas(n));
    load(fullfile(outDir,outputDir,dataFile))
    D = outputs.D;
    N = outputs.N;
    M = outputs.M;
    y = outputs.y;
    X = outputs.X;
    scales = outputs.scales;
    center = (M+1)/2;
    
    AD = reSampleCustomArrayCenter(N,D,scales,center);
    AD = padarray(AD,[0 M-1 0],0,'post');
    ADf = fft2(AD);
    Yhat = unpad(squeeze(ifft2(sum(bsxfun(@times,ADf,fft2(X)),3),'symmetric')),M-1,'pre');
    Yhat = gather(Yhat);

    trueErr(n) = norm(y_true-Yhat);
    dataErr(n) = norm(squeeze(outputs.y)-Yhat);
    noiseNorm(n) = norm(randn(N,T)*sigmas(n));
    noiseNorm2(n) = norm(y_true-squeeze(outputs.y));

    ff = figure();
    hold on
    kk = 1;
    for ttt = [1,11,21]
        subplot(3,1,kk)
        plot(outputs.y(:,:,ttt),'-')
        hold on
        plot(y_true(:,ttt),'-')
        plot(Yhat(:,ttt),'-')
        kk = kk + 1; 
    end
    ff.Position = [941 59 560 825];
    legend('data','truth','recon')
    pause()
end

%%
figure()
hold on
plot(meanSNR(2:NN),trueErr(2:NN),'o-')
plot(meanSNR(2:NN),dataErr(2:NN),'x-')
plot(meanSNR(2:NN),noiseNorm(2:NN),'s-')
plot(meanSNR(2:NN),noiseNorm2(2:NN),'x-')
xlabel('SNR','Fontsize',14)
ylabel('Error','Fontsize',14)
legend('$\|\hat{{\bf b}}-{\bf f}\|_2$','$\|\hat{{\bf b}}-{\bf b}\|_2$',...
    '$\|{\bf w}\|_2$','$\|{\bf b}-{\bf f}\|_2$',...
    'interpreter','latex','Fontsize',14)
%% Next copy figures associated with selected parameters to a folder

% for n = 2:NN
%     outputDir = "gaus_example_8_23_24_X0_D0_V00_sig_"+num2str(n);
% end
