%% Load image

% Data I/O directories
data_dir = 'E:\CHESS_data\al7075_311_polar\';
results_dir = 'E:\CHESS_results\test\';

% data_dir = 'E:\MMPAD_data\ring1_zero\';
% results_dir = 'E:\MMPAD_data\test_output\';

% Select image
P.set = 1;
P.img = 50;

% Load polar_image
load([data_dir,... 
'polar_image_',...
num2str(P.set),'_',...
num2str(P.img),'.mat']);

B = polar_image(:,600:900);

% Load image prefit to use as weights
% load('.mat');

% Ring sampling parameters
P.num_theta= size(B,2);
P.num_rad = size(B,1);
P.sampleDims = [546,1];

% Basis function variance parameters
P.basis = 'norm2';
P.num_var_t = 20;
P.num_var_r = 8;
P.var_theta = linspace(0.5,  20,P.num_var_t).^2;
P.var_rad = linspace(0.5,      6,P.num_var_r).^2;

% Zero padding and mask
% maskCols = 129:133;
% zPad = [0,0];
% zMask = zeros(size(zeroPad(polar_image,zPad)));
% zMask(:,maskCols) = 1;
% zMask = onePad(zMask,zPad);
% [r,c] = find(zMask==1);
% zMask = [r,c];

zPad = [0,0];
zMask = [];

% Create dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack(P);

%Show dictionary atoms
% P_plot = P;
% P_plot.basis = 'max';
% A0_stack = unshifted_basis_matrix_stack(P_plot);
% uplim = max(A0_stack(:));
% 
% dict_fig = figure(2);
% [ha, pos] = tight_subplot(P.num_var_t/2,P.num_var_r,[0.01 0.01],[.02 .08],[.01 .01]); 
% k = 1;
% for i = 2:2:P.num_var_t
%     for j = 1:P.num_var_r
%         axes(ha(k))
%         center = shift2D(A0_stack(:,:,i,j),floor(P.num_rad/2),floor(P.num_theta/2));
%         imshow(center(:,100:200),'DisplayRange',[0 uplim],'Colormap',jet)
%         k = k + 1;
%     end
% end

x = B;
params.epsilon = 0.2;
params.delta = 5e-4;
params.showImage = 1;
params.isNonnegative = 1;
params.zeroMask = zMask;
%% Group Convolutional Matching Pursuit (GCMP)
a = CMP(A0ft_stack,x,params);

%% Compare test images in terms of AWMV

% GCMP solutions
% Ring sampling parameters
% P.num_theta= size(polar_image,2);
% P.num_rad = size(polar_image,1);
% P.dtheta = 1;
% P.drad = 1;
% 
% % Basis function dictionary parameters
% P.num_var_t = 12;
% P.num_var_r = 8;
% P.var_theta = linspace(1, 20, P.num_var_t).^2;
% P.var_rad   = linspace(1, 4, P.num_var_r).^2;
% A0ft_stack = unshifted_basis_matrix_ft_stack(P);
% for num = [0,4]
%     load(sprintf('al7075_gcmp_%i_35.mat',num))
%     sparsity = sum(a(:)>0);
%     fit = Ax_ft_2D(A0ft_stack,a);
%     err = norm(fit(:)-polar_image(:))/norm(polar_image(:));
%     awmv = computeAWMV(a,sqrt(P.var_theta),sqrt(P.var_rad));
%     fprintf(['At load %i... \n',...
%              'sparsity = %i \n',...
%              'Relative error = %0.3f \n',... 
%              'AWMV = %2.3f \n\n'],num,sparsity,err,awmv)
%          
% figure(num+1)
% subplot(3,1,1)
% imshow(polar_image,'DisplayRange',[0 900],'Colormap',jet)
% subplot(3,1,2)
% imshow(fit,'DisplayRange',[0 900],'Colormap',jet)
% subplot(3,1,3)
% imshow(abs(polar_image-fit),'DisplayRange',[0 900],'Colormap',jet)
% end


