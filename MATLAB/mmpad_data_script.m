%% Load data into single array
img_array = zeros(67,265,396);
for i = 1:67
    fname = sprintf('mmpad_img_%i.mat',i-1);
    load(fullfile('E:\PureTiRD_nr2_c_x39858',fname))
    img_array(i,:,:) = squeeze(sum(mmpad_img,1));
end

%% Crop image array to separate rings
% rows 1,130-134
% These remove the zero pixels
% ring1 = img_array(:,[2:129,135:262],5:40);
% ring2 = img_array(:,[2:129,135:262],210:260);
% ring3 = img_array(:,[2:129,135:262],268:310);
% ring4 = img_array(:,[2:129,135:262],355:395);
ring1 = img_array(:,2:262,5:40);
ring2 = img_array(:,2:262,210:260);
ring3 = img_array(:,2:262,268:310);
ring4 = img_array(:,2:262,355:395);
% Set ring
ring = ring3;

%% Save reduced images to file
imDir = 'D:\MMPAD_data\';
for ring_num = 3
    fdir = sprintf('ring%i_zero',ring_num);
    mkdir(fullfile(imDir,fdir))
    for i = 1:546
        fname = sprintf('mmpad_img_%i.mat',i);
        polar_image = squeeze(ring(i,:,:));
        save(fullfile(imDir,fdir,fname),'polar_image')
    end
end


%% Display images
for j = 1:1:546
    figure(1)
    subplot(1,3,1)
    imshow(squeeze(ring(j,:,:)),'DisplayRange',[0,10e3],'ColorMap',jet)
    subplot(1,3,2)
    imshow(squeeze(omega(j,:,:)),'DisplayRange',[0,1],'ColorMap',jet)
%     subplot(1,3,3)
%     imshow(squeeze(ring_fill(j,:,:)),'DisplayRange',[0,10e3],'ColorMap',jet)
     pause
end

%% Load saved array
load('D:\MMPAD_processing\mmpad_summed.mat')


%% Display array
for j = 1:546
    imshow(squeeze(ring_fill(j,:,:)),'DisplayRange',[0,10e3],'ColorMap',jet)
    pause(0.1)
end

%% Tensor completion on array

% Setup data
normalize = max(abs(ring(:)));
Tn = ring/normalize;
A = diag(sparse(double(omega(:))));
b = A*Tn(:);
[n1,n2,n3] = size(Tn);

% Parameters
rho = 10;
alpha = 1;
myNorm = 'tSVD_1';
maxIter = 1000;
[X,history]  =    tensor_cpl_admm( A , b , rho , alpha , ...
                     [n1,n2,n3] , maxIter , myNorm , 0 );  

%% Display new array
ring_fill = reshape(X,[n1,n2,n3])*normalize; 
%% Save new array
save('D:\MMPAD_processing\mmpad_summed.mat','img_array')

%% Test FISTA 
load('D:\MMPAD_data\ring1_zero\mmpad_img_99.mat')
%polar_image = polar_image(224:end,:);

% Zero pad and Zero mask
% ring1 = img_array(:,[2:129,135:262],5:40);
% ring2 = img_array(:,[2:129,135:262],210:260);
% ring3 = img_array(:,[2:129,135:262],268:310);
% ring4 = img_array(:,[2:129,135:262],355:395);

% Zero padding and mask
maskCols = 129:133;
zPad = [0,0];
zMask = zeros(size(zeroPad(polar_image,zPad)));
zMask(:,maskCols) = 1;
zMask = onePad(zMask,zPad);
[r,c] = find(zMask==1);
zMask = [r,c];

% Ring sampling parameters
P.num_theta= size(polar_image,2);
P.num_rad = size(polar_image,1);
P.dtheta = 1;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 8;
P.num_var_r = 6;
P.var_theta = linspace(P.dtheta/2,32, P.num_var_t).^2;
P.var_rad   = linspace(P.drad/2,  6, P.num_var_r).^2;

% Generate unshifted basis function matrices
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
A0_stack = unshifted_basis_matrix_stack_norm2(P);

% FISTA parameters
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 500;
params.lambda = 0.01;
params.gamma = 1;
params.beta = 1.1;
params.maxIter = 500;
params.maxIterReg = 500;
params.isNonnegative = 1;
params.zeroPad = zPad;
params.zeroMask = zMask;
params.noBacktrack = 0;
params.plotProgress = 0;
P.params = params;

% Initialize solution
x_init = rand(size(A0ft_stack));
x_init = x_init/norm(x_init(:));

% FISTA with backtracking
polar_image = polar_image/norm(polar_image(:));
[x_hat, err, obj, l_0] = FISTA_Circulant(A0ft_stack,polar_image,x_init,params);

%  save(sprintf('fista_fit_%i_%i.mat','x_hat','err','polar_image','P',...
%       P.set,P.img))
%% View mmpad fit
img_fit = Ax_ft_2D(A0ft_stack,x_hat);
img_fit_z = forceMaskToZero(img_fit,zMask);
rel_err = norm(polar_image(:)-img_fit_z(:))/norm(polar_image(:))
l0 = sum(x_hat(:)>0)
lim1 = 0;
lim2 = max(polar_image(:));
% Plot both images
figure(6)
subplot(1,3,1)
imshow(polar_image,'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,3,2)
imshow(img_fit,'DisplayRange',[lim1 lim2],'Colormap',jet);
subplot(1,3,3)
imshow(abs(polar_image-img_fit),'DisplayRange',[lim1 2*lim2],'Colormap',jet);

%% View basis functions used
figure(7)
var_signal = squeeze(sum(sum(x_hat,1),2));
lim2 = max(var_signal(:));
imshow(var_signal,'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');

%%
lim2 = max(x_hat(:));
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        subplot(1,P.num_var_r,j)
        imshow(x_hat(:,:,i,j),'DisplayRange',[lim1 lim2],'Colormap',jet,'InitialMagnification','fit');
    end
    pause
end

        %% View basis functions
figure(5)
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        subplot(P.num_var_t,P.num_var_r,(i-1)*P.num_var_r+j)
        basis_func = shift2D(squeeze(A0_stack(:,:,i,j)),round(P.num_rad/2),round(P.num_theta/2));
        imshow(basis_func','DisplayRange',[lim1 0.8],'Colormap',jet);
    end
end

%% Test basis functions
x_test = zeros(size(A0ft_stack));
x_test(15:7:end,1:5:10,1,1)=1;
test_img = Ax_ft_2D(A0ft_stack,x_test);
lim1 = 0;
lim2 = max(test_img(:));
figure(8)
imshow(test_img,'DisplayRange',[lim1 lim2],'Colormap',jet);

