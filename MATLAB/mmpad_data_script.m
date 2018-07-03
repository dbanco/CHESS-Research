%% Load data into single array
% img_array = zeros(546,265,396);
% for i = 1:546
%     fname = sprintf('mmpad_img_%i.mat',i-1);
%     load(fullfile('D:\MMPAD_data',fname))
%     img_array(i,:,:) = squeeze(sum(mmpad_img,1));
% end

%% Crop image array to separate rings
% rows 1,130-134
ring1 = img_array(:,[2:129,135:262],5:40);
ring2 = img_array(:,[2:129,135:262],210:260);
ring3 = img_array(:,[2:129,135:262],268:310);
ring4 = img_array(:,[2:129,135:262],355:395);

% Set ring
ring = ring2;

omega = ones(size(ring));
omega(:,130:134,:) = 0;

%% Display iamges
for j = 1:1:546
    figure(1)
    subplot(1,3,1)
    imshow(squeeze(ring(j,:,:)),'DisplayRange',[0,10e3],'ColorMap',jet)
    subplot(1,3,2)
    imshow(squeeze(omega(j,:,:)),'DisplayRange',[0,1],'ColorMap',jet)
%     subplot(1,3,3)
%     imshow(squeeze(ring_fill(j,:,:)),'DisplayRange',[0,10e3],'ColorMap',jet)
     pause(0.1)
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

%% Test FISTA Masked
% Ring sampling parameters
P.ring_width = 20;
P.num_theta= 396;
P.num_rad = 265;
P.dtheta = 1;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 10;
P.num_var_r = 10;
P.var_theta = linspace(P.dtheta,30,  P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  100, P.num_var_r).^2;

% Generate unshifted basis function matrices
P.betap = P.dtheta*P.drad;
P.weight = 1;
P.alphap = 10;

A0ft_stack = unshifted_basis_matrix_ft_stack(P);
A0_stack = unshifted_basis_matrix_stack(P);

% FISTA parameters
params.stoppingCriterion = 1;
params.tolerance = 1e-6;
params.L = 1e1;
params.lambda = 50;
params.beta = 1.2;
params.maxIter = 100;
params.isNonnegative = 1;
P.params = params;

% Set image
test_im = squeeze(img_array(1,:,:));

[i1,i2]=find(test_im == 0);
mask = [i1,i2];

x_init = ones(size(A0ft_stack));
%% FISTA with backtracking
[x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,test_im,x_init,params);

err(end)
%  save(sprintf('fista_fit_%i_%i.mat','x_hat','err','polar_image','P',...
%       P.load_step,P.img))