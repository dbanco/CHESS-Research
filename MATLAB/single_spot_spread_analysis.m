%% Load data
data1 = load(fullfile('D:','Chess_data','al7075_222_polar_fit2','fista_fit_0_175_task_175.mat'));
data2 = load(fullfile('D:','Chess_data','al7075_222_polar_fit2','fista_fit_4_175_task_995.mat'));
A0ft_stack = unshifted_basis_matrix_ft_stack(P);

y_fit1 = Ax_ft_2D(A0ft_stack,data1.x_hat);
y_fit2 = Ax_ft_2D(A0ft_stack,data2.x_hat);

%% Script to develop spot-by-spot analysis
figure(2)
subplot(3,1,1)
imshow(real(log(data1.polar_image)),'DisplayRange',[0 9],'Colormap',jet)
title('Data')

subplot(3,1,2)
imshow(real(log(data2.polar_image)),'DisplayRange',[0 9],'Colormap',jet)
title('Data')

subplot(3,1,3)
imshow(real(log(y_fit1)),'DisplayRange',[0 9],'Colormap',jet)
title('Fit')

% Sum coefficients (sum over time as well)
x_hat_sum1 = squeeze(sum(sum(data1.x_hat,3),4));
x_hat_sum2 = squeeze(sum(sum(data2.x_hat,3),4));
x_hat_sum = x_hat_sum1 + x_hat_sum2;
% Convert to binary
x_hat_bin = x_hat_sum > 0;

% plot binary
maxLim = log(max(x_hat_sum(:)));
subplot(3,1,3)
imshow(x_hat_bin)
%imshow(real(log(x_hat_sum)),'DisplayRange',[0 maxLim],'Colormap',jet)
title('Binary')

%% Compute metric for each single spot
% Get connected components
x_cc_struct = bwconncomp(x_hat_bin);
vdf1 = {};
vdf2 = {};
for i = 1:numel(x_cc_struct.PixelIdxList)
    % Get spot
    pixels = x_cc_struct.PixelIdxList{i};
    % Convert to 2d pixel indices
    [j,k] = ind2sub(size(polar_image),pixels);
    % Compute variance distribution functions at both times
    vdf1{i} = squeeze(sum(sum(data1.x_hat(j,k,:,:),1),2));
    vdf2{i} = squeeze(sum(sum(data2.x_hat(j,k,:,:),1),2));
    vdf1{i} = vdf1{i}./sum(vdf1{i}(:));
    vdf2{i} = vdf2{i}./sum(vdf2{i}(:));
end

%% View vdfs
n = 9;
maxLim = max(max(vdf1{n}(:)),max(vdf2{n}(:)));
figure(3)
subplot(1,2,1)
imshow(vdf1{n},'DisplayRange',[0 maxLim],'Colormap',jet)
title('Initial')

subplot(1,2,2)
imshow(vdf2{n},'DisplayRange',[0 maxLim],'Colormap',jet)
title('Loaded')

