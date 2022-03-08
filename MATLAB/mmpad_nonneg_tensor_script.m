%% Load MMPAD sequence
r = 1;
data_dir = ['D:\MMPAD_data_nr1\ring', num2str(r), '_zero'];
T = 200;

for t = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t),'.mat']))
    if t == 1
        [N,M] = size(polar_image);
        X = zeros(N,M,T);
    end
    X(:,:,t) = polar_image;
end
Xt = tensor(X);

%% Nonnegative CP with Dictionary learning

% Initialize NCP decomp
R = 20;
F = ncp(Xt,R,'method','hals');

% Initialize dictionaries
P.params.zeroPad = [];
P.basis = 'norm2';
K = 15;
P.num_var_t = K;
P.num_theta = N;
P.var_theta = linspace(0.5,5,P.num_var_t).^2;
A0ft_rad = unshifted_basis_vector_ft_stack_zpad(P);

P2 = P;
P2.num_theta = M;
P2.var_theta = linspace(0.5,50,P.num_var_t).^2;
A0ft_az = unshifted_basis_vector_ft_stack_zpad(P2); 

% Initialize sparse coding
params.lambda1 = 0.2;
params.rho1 = 1;
params.adaptRho = 1;
params.tau = 1.05;
params.mu = 2;
params.alpha = 1.8;
params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 800;
params.tolerance = 1e-8;
params.isNonnegative = 1;
params.zeroPad = [];
params.zeroMask = [];
params.plotProgress = 0;
params.verbose = 0;
X1 = zeros(N,K,R);
X2 = zeros(M,K,R);
Y = F.U;
for r = 1:R
    sprintf('%i of %i\n',r,R)
    X1(:,:,r) = convADMM_LASSO_Sherman_1D(A0ft_rad,F.U{1}(:,r),zeros(N,15),params);
    X2(:,:,r) = convADMM_LASSO_Sherman_1D(A0ft_az, F.U{2}(:,r),zeros(M,15),params);
    Y{1}(:,r) = Ax_ft_1D(A0ft_rad,X1(:,r));
    Y{2}(:,r) = Ax_ft_1D(A0ft_az,X2(:,r));
end

%% Run constrained HALS
F1 = ncp(Xt,R,'method','hals_dict','stepSize',1,'Y',Y);


%% Recon
F_recon = double(tensor(F1));

%% View result
figure(1)
n = 20;

subplot(2,1,1)
imagesc(F_recon(:,:,n))
title('recon')

subplot(2,1,2)
imagesc(X(:,:,n))
title('data')

%% Visualize factors
U1 = F.U{1};
U2 = F.U{2};
U3 = F.U{3};
% 
% figure(2) % View angle components
% for i = 1:20
%     subplot(5,4,i)
%     angComp = F.lambda(i)*U1(:,i)*U2(:,i)';
%     imagesc(angComp)
% end
% 
% figure(3) % View radial/time components
%    for i = 1:20
%     subplot(5,4,i)
%     radComp = F.lambda(i)*U1(:,i)*U3(:,i)';
%     imagesc(radComp)
% end
% 
% figure(4) % View az/time components
%    for i = 1:20
%     subplot(5,4,i)
%     azComp = F.lambda(i)*U2(:,i)*U3(:,i)';
%     imagesc(azComp)
%    end

%% Sort by time occurence
% Compute time center of mass
timeW = 1:T;
timeCoM = zeros(20,1);
for i = 1:20
    timeCoM(i) = timeW*U3(:,i)/sum(U3(:,i),'all');
end
[~,ind] = sort(timeCoM);

%% Visualize factors in time order
U1 = F.U{1};
U2 = F.U{2};
U3 = F.U{3};

figure(2) % View angle components
for i = 1:20
    subplot(5,4,i)
    j = ind(i);
    angComp = F.lambda(j)*U1(:,j)*U2(:,j)';
    imagesc(angComp)
    title(['t= ',num2str(timeCoM(j))])
end

figure(3)
U3sort = U3;
for i = 1:20
    subplot(5,4,i)
    j = ind(i);
    plot(U3(:,j))
    U3sort(:,i) = U3(:,j);
end

figure(4)
imagesc((U3sort'))

figure(5)
for i = 1:20
    subplot(5,4,i)
    j = ind(i);
    k = floor(timeCoM(j));
    imagesc(X(:,:,k))
    title(['t= ',num2str(k)])
end

%% 2D CSC for R components
angComp = U1(:,1)*U2(:,1)';
%Parameters
P.num_rad = size(angComp,1);
P.num_theta = size(angComp,2);

% Basis function variance parameters
P.num_var_t = 15;
P.num_var_r = 10;
P.var_theta = linspace(1,50,P.num_var_t).^2;
P.var_rad   = linspace(1, 3,    P.num_var_r).^2;
P.basis = 'norm2';

% fista params
params.lambda1 = 0.01; % sparsity penalty
params.rho1 = 1;  % initial ADMM
params.adaptRho = 1; % binary flag for adaptive rho
params.mu = 10;       % tolerated factor between primal/dual residual
params.tau = 1.05;   % rho update factor
params.alpha = 1.8; % over-relaxation paramter
params.isNonnegative = 1; % flag to enforce nonnegativity
params.stoppingCriterion = 'OBJECTIVE_VALUE';
params.maxIter = 100;
params.tolerance = 1e-6;
zPad = [0,0];
zMask = [];
params.zeroPad = zPad; % number of [rows,columns]of padding to add to data
params.zeroMask = zMask; % specifies columns known to be zero

params.plotProgress = 0; % flag to plot intermediate solution at each iteration 
params.verbose = 1;      % flag to print objective values at each iteration 
P.params = params;

A0ft_stack = unshifted_basis_matrix_ft_stack(P);

% Initialize solution
x_init = zeros(size(A0ft_stack));
b = polar_image/100;
[n1,n2] = size(b);
% Solve
% [x_hat,err,obj] = convADMM_LASSO_Sherman_2D(A0ft_stack,b,x_init,params);
params.L = 1;
params.lambda = 0.001;
params.beta = 2;
params.noBacktrack = 0;
R = 20;
X_hat = zeros(n1,n2,P.num_var_t,P.num_var_r,R);
for i = 1:R
    j = ind(i);
    angComp = U1(:,j)*U2(:,j)';
    [x_hat,err,obj] = FISTA_Circulant(A0ft_stack,angComp,x_init,params);
    X_hat(:,:,:,:,i) = x_hat;
end

%% Compute AWMV for R components
awmv_az = zeros(1,R);
awmv_rad = zeros(1,R);
U3sort = U3;
for i = 1:R
	[az,rad] = computeAWMV(X_hat(:,:,:,:,i),...
               sqrt(P.var_theta),sqrt(P.var_rad) );
    awmv_az(i) = az;
    awmv_rad(i) = rad;
    j = ind(i);
    U3sort(:,i) = U3(:,j);
end
U3scale = U3sort;
for t = 1:T
   U3scale(t,:) = U3sort(t,:)/sum(U3sort(t,:));
end
figure(6)
imagesc(U3scale)
full_awmv_az = awmv_az*U3scale';
plot(full_awmv_az)
hold on 
for i = 1:20
    j = ind(i);
    k = floor(timeCoM(j));
   stem(k,awmv_az(i),'o','MarkerSize',4)
end
title('AWMV computed from Nonneg CP')
xlabel('time')
ylabel('AWMV')

%% Nonneg CP Error and CSC Error
rel_cp_err = zeros(1,T);
rel_csc_err = zeros(1,T);
U1xU2 = zeros(n1,n2,R);
for i = 1:R
    U1xU2(:,:,i) = F.lambda(i)*Ax_ft_2D(A0ft_stack,X_hat(:,:,:,:,i));
    j = ind(i);
    U3sort(:,i) = U3(:,j);
end
F_csc = double(ttm(tensor(U1xU2),U3sort,3));
for t = 1:T
%     F_csc(:,:,t) = U1xU2*U3sort(t,:)';
    rel_cp_err(t) = norm(X(:,:,t)-F_recon(:,:,t))/norm(X(:,:,t));
    rel_csc_err(t) = norm(X(:,:,t)-F_csc(:,:,t))/norm(X(:,:,t));
end

figure(7)
plot(rel_cp_err)
hold on
plot(rel_csc_err)
legend('Error due to CP','Error after CSC')
ylabel('Relative Error')
xlabel('time')

t=500;
figure(8)
subplot(3,1,1)
imagesc(X(:,:,t))
title(['data t = ',num2str(t)])
subplot(3,1,2)
imagesc(F_recon(:,:,t))
title('CP Recon')
subplot(3,1,3)
imagesc(F_csc(:,:,t))
title('CP+CSC Recon')
