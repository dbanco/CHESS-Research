%% SVD Analysis of BCCB matrix
%% Fixed Parameters
clear all; close all
% Ring sampling parameters
P.ring_width = 2;
P.num_theta = 11;
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

A0_stack = unshifted_basis_matrix_stack(P);
basis1 = squeeze(A0_stack(:,:,1,1));
basis2 = squeeze(A0_stack(:,:,15,10));

basis = basis2;

[n,m] = size(basis); 
basis_shift_row = [basis(n+1-3:n,:); basis(1:n-3,:)];
% column shift
basis = [basis_shift_row(:,m+1-6:m), basis_shift_row(:,1:m-6)];

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

bccb2 = bccb;

figure(2)
imshow(bccb,'DisplayRange',[0 1],'ColorMap',jet)
colorbar()
title('BCCB')
