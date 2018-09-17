%% Fixed Parameters
clear all; close all
% Ring sampling parameters
P.ring_width = 3;
P.num_theta = 18;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 2*pi/P.num_theta;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;

P.var_theta = linspace(P.dtheta,pi/3,P.num_var_t).^2;
P.var_rad   = linspace(P.drad,  1.3,    P.num_var_r).^2;

% basis weighting
P.weight = 0;
P.betap = P.dtheta*P.drad;
P.alphap = 1;

A0ft_stack = unshifted_basis_matrix_ft_stack(P);
A0_stack = unshifted_basis_matrix_stack(P);
basis = squeeze(A0_stack(:,:,1,1));
[n,m] = size(basis); 
% basis_shift_row = [basis(n+1-3:n,:); basis(1:n-3,:)];
% % column shift
% basis = [basis_shift_row(:,m+1-6:m), basis_shift_row(:,1:m-6)];

% [n,m] = size(basis1); 
% basis = [flipud(fliplr(basis1(2:n,2:m))),flipud(basis1(2:n,1:m));
%          fliplr(basis1(1:n,2:m)), basis1];
% basis = [basis1,zeros(n,1), fliplr(basis1(1:n,1:m-1));
%          zeros(1,2*m-1);
%          flipud(basis1(1:n-1,1:m)),zeros(n-1,1), flipud(fliplr(basis1(1:n-1,1:m-1)))];

figure(1)
imshow(basis,'DisplayRange',[0 1],'ColorMap',jet)

% Concatenate 2D shifted basis functions to form BCCB matrix
[n,m] = size(basis); 
N = numel(basis(:));
bccb = zeros(N,N);
k=1;
for j = 0:size(basis,2)-1
    for i = 0:size(basis,1)-1
        % row shift
        basis_shift_row = [basis(n+1-i:n,:); basis(1:n-i,:)];
        % column shift
        basis_shift_row_col = [basis_shift_row(:,m+1-j:m), basis_shift_row(:,1:m-j)];
        % vectorize
        basis_vec = basis_shift_row_col(:);
        % add to matrix
        bccb(:,k) = basis_vec;
        k = k+1;
    end
end
% figure(2)
% imshow(bccb,'DisplayRange',[0 1],'ColorMap',jet)
% colorbar()
% title('BCCB')

% figure(3)
% for k = 1:m
%     subplot(3,4,k)
%     imshow(bccb(n*(k-1)+1:n*k,1:n),'DisplayRange',[0 1],'ColorMap',jet)
%     title(sprintf('C_{%i}',k-1))
% end

%% Test bccb matrix
x = rand(P.num_rad,P.num_theta,1,1);
x_vec = squeeze(x(:));
y1 = Ax_ft_2D(A0ft_stack(:,:,1,1),x);
y2 = bccb*x_vec;

err1 = norm(y1(:)-y2)

figure(2)
subplot(2,1,1)
imshow(reshape(y2,[P.num_rad,P.num_theta]),'DisplayRange',[min(y2(:)) max(y2(:))],'ColorMap',jet)
colorbar()

subplot(2,1,2)
imshow(reshape(y1,[P.num_rad,P.num_theta]),'DisplayRange',[min(y1(:)) max(y1(:))],'ColorMap',jet)
colorbar()

r = rand(P.num_rad,P.num_theta);
r_vec = squeeze(r(:));
g1 = bccb'*r_vec;
g2 = AtR_ft_2D((A0ft_stack(:,:,1,1)),r);

% g2 = shift2D(g2,2,2);

err2 = norm(g1(:)-g2(:))

figure(3)
subplot(2,1,1)
imshow(reshape(g2(:,:,1,1),[P.num_rad,P.num_theta]),'DisplayRange',[min(g2(:)) max(g2(:))],'ColorMap',jet)
colorbar()

subplot(2,1,2)
imshow(reshape(g1(:,:,1,1),[P.num_rad,P.num_theta]),'DisplayRange',[min(g1(:)) max(g1(:))],'ColorMap',jet)
colorbar()

figure(4)
subplot(2,1,1)
imshow(bccb,'DisplayRange',[min(bccb(:)) max(bccb(:))],'ColorMap',jet)
colorbar()

bccbt = bccb';
subplot(2,1,2)
imshow(bccbt,'DisplayRange',[min(bccb(:)) max(bccb(:))],'ColorMap',jet)
colorbar()

figure(5)
subplot(2,1,1)
t1 = 17;
imshow(reshape(bccb(:,t1),[P.num_rad,P.num_theta]),'DisplayRange',[min(bccb(:,t1)) max(bccb(:,t1))],'ColorMap',jet)
colorbar()
bccbt = bccb';
subplot(2,1,2)
t2 = 1;
imshow(reshape(bccbt(:,t2),[P.num_rad,P.num_theta]),'DisplayRange',[min(bccb(:,t2)) max(bccb(:,t2))],'ColorMap',jet)
colorbar()
