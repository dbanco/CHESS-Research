
clear all
close all


%% Define parameters

% Length of intensity data (theta coordinate)
P.num_theta = 30; 

% Define dictionary of Gaussian basis functions
P.num_var_t = 1;   % Number of different basis functions 
P.var_theta = linspace(1/2,50,P.num_var_t).^2; % Variances of basis functions
P.basis = 'norm2';
P.params.zeroPad = 0;

K = P.num_var_t;
N = P.num_theta;

% Construct dictionary
A = zeros(P.num_theta,P.num_theta*P.num_var_t);
A0ft_stack = unshifted_basis_vector_ft_stack_zpad(P);
A0_stack = unshifted_basis_vector_stack_zpad(P);
for j = 1:K
    for i = 1:P.num_theta
        ind1 = 1 + (i-1)*N + (j-1)*K*N;
        ind2 =         i*N + (j-1)*K*N;
        ind3 = i + (j-1)*N;
        A(:,ind3) = shift1D(A0_stack(:,j),i-1);
    end
end
x = randn(P.num_theta,P.num_var_t);

y1 = Ax_ft_1D(A0ft_stack,x);
y2 = A*x(:);

norm(y1-y2)/norm(y2)

Aty1 = AtR_ft_1D(A0ft_stack,y1);
Aty2 = A'*y1;

norm(Aty1(:)-Aty2)/norm(Aty2)

figure(1)
subplot(2,1,1)
hold on
plot(y1,'o')
plot(y2)

subplot(2,1,2)
hold on
plot(Aty1(:),'o')
plot(Aty2(:))