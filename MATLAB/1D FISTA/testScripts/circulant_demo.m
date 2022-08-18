%% Fixed Parameters
clear all; close all
% Ring sampling parameters
P.ring_width = 3;
P.num_theta = 18;
P.num_rad = 2*P.ring_width+1;
P.dtheta = 1;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 15;

P.var_theta = linspace(P.dtheta,pi/3,P.num_var_t).^2;

% basis weighting
P.weight = 0;
P.betap = P.dtheta*P.drad;
P.alphap = 1;

A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
A0_stack = unshifted_basis_vector_stack_norm2(P);
basis = squeeze(A0_stack(:,1));
[n] = numel(basis); 
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
plot(basis)

circulant = zeros(n,n);
k=1;

for i = 0:size(basis,1)-1
    % row shift
    basis_shift_row = shift1D(basis,i);
    % add to matrix
    circulant(:,k) = basis_shift_row';
    k = k+1;
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
x = rand(P.num_theta,1);
x_vec = squeeze(x(:));
y1 = Ax_ft_1D(A0ft_stack(:,1),x);
y2 = circulant*x_vec;

err1 = norm(y1(:)-y2)

figure(2)
subplot(2,1,1)
plot(y2)

subplot(2,1,2)
plot(y1)

r = rand(P.num_theta,1);
r_vec = squeeze(r(:));
g1 = circulant'*r_vec;
g2 = AtR_ft_1D((A0ft_stack(:,1)),r);

% g2 = shift2D(g2,2,2);

err2 = norm(g1(:)-g2(:))

figure(3)
subplot(2,1,1)
plot(g2)

subplot(2,1,2)
plot(g1)

figure(4)
subplot(2,1,1)
imshow(circulant,'DisplayRange',[min(circulant(:)) max(circulant(:))],'ColorMap',jet)
colorbar()

bccbt = circulant';
subplot(2,1,2)
imshow(bccbt,'DisplayRange',[min(circulant(:)) max(circulant(:))],'ColorMap',jet)
colorbar()

figure(5)
subplot(2,1,1)
t1 = 17;
plot(circulant(:,t1))
colorbar()
bccbt = circulant';
subplot(2,1,2)
t2 = 1;
plot(bccbt(:,t2))
colorbar()
