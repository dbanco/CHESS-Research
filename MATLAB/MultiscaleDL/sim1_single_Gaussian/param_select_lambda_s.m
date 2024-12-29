function [lambda_s,outInd] = param_select_lambda_s(outputDir,tradeoff,scaleP,fig_num,criterion,sigma,y_true)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

files = dir(fullfile(outputDir, '*.mat'));

% Extract the file names and store them in a cell array
matFileNames = {files.name};

NN = numel(matFileNames);
rel_error = zeros(NN,1);
true_error = zeros(NN,1);
l1_norm = zeros(NN,1);
lambda_s_vec = zeros(NN,1);
for i = 1:numel(matFileNames)
    % Load outputs
    load(fullfile(outputDir,matFileNames{i}))

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
    true_error(i) = sqrt(sum((y_true-Yhat).^2,'all'));
    % Compute L1-norm
    l1_norm(i) = norm(X(:),1);
    % Get lambda parameter
    lambda_s_vec(i) = outputs.lambda;
end

[lambda_s_sort,ind] = sort(lambda_s_vec);
true_sort = true_error(ind);
err_sort = rel_error(ind);
l1_sort = l1_norm(ind);

% Normalized origin distance criterion
switch criterion
    case 'origin_dist'
        crit= tradeoff*abs((err_sort-scaleP(1))/scaleP(2)) +...
                        abs((l1_sort-scaleP(3))/scaleP(4));
        [~, selInd] = min(crit);
        lambda_s = lambda_s_sort(selInd);
    case 'discrepancy'
        crit = abs(err_sort/sqrt(N*T) - sigma);
        [~,selInd] = min(crit);
        lambda_s = lambda_s_sort(selInd);
    case 'truth_error'
        [~,selInd] = min(true_sort);
        lambda_s = lambda_s_sort(selInd);
    case 'curvature'
        % Curvature criterion
        [b, a] = butter(4, 0.1, 'low');  % low-pass Butterworth filter
        z = filtfilt(b,a,log(err_sort));
        x = filtfilt(b,a,log(l1_sort));
        
        dlam = diff(lambda_s_sort);
        dz = diff(z)./dlam;
        dx = diff(x)./dlam;
        ddz = diff(z,2)./dlam(1:end-1);
        ddx = diff(x,2)./dlam(1:end-1);
        
        curvature = (dz(1:end-1).*ddx - dx(1:end-1).*ddz)./...
                    (dx(1:end-1).^2 + dz(1:end-1).^2).^1.5;
        
        [~, selInd] = max(curvature);
        lambda_s = lambda_s_sort(selInd);
    case 'curvature_poly'
        % Curvature criterion
        z = log(err_sort);
        x = log(l1_sort);
        control_x = polyfilt(x,lambda_s_sort);
        control_z = polyfilt(z,lambda_s_sort);
        control_lam = lambda_s_sort(3:end-2);
        xq = linspace(min(x),max(x),100);
        zq = spline(control_x,control_z,xq);
%         zq2 = spline(x,z,xq);
        figure
        hold on
        plot(x,z,'o-',control_x,control_z,'x-',xq,zq)
        
%         dlam = diff(control_lam);
%         dz = diff(zq2)./dlam;
%         dx = diff(xq)./dlam;
%         ddz = diff(zq2,2)./dlam(1:end-1);
%         ddx = diff(xq,2)./dlam(1:end-1);

        dz = diff(zq);
        dx = diff(xq);
        ddz = diff(zq,2);
        ddx = diff(xq,2);
        
        curvature = (dz(1:end-1).*ddx - dx(1:end-1).*ddz)./...
                    (dx(1:end-1).^2 + dz(1:end-1).^2).^1.5;
        figure
        plot(curvature)
        [~, selInd] = max(curvature);
        lambda_s = lambda_s_sort(selInd);
end
if nargin > 3
    figure(fig_num)
    loglog(l1_sort(1:end),err_sort(1:end),'o-')
    ylabel('Error')
    xlabel('l_1-norm')
    hold on
    loglog(l1_sort(selInd),err_sort(selInd),'sr','MarkerSize',10)
end

outInd = ind(selInd);

end