%% Define spots
% Ring sampling parameters
ring_width = 15;
P.num_theta= 60;
P.num_rad = 2*ring_width+1;
P.dtheta = 1;
P.drad = 1;

% Basis function variance parameters
P.num_var_t = 4;
P.num_var_r = 3;
P.var_theta = linspace(P.dtheta, 16, P.num_var_t).^2;
P.var_rad   = linspace(P.drad,   3, P.num_var_r).^2;

% Create dictionary
A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(P);
A0_stack = unshifted_basis_matrix_stack_norm2(P);
norms = zeros(P.num_var_t,P.num_var_r);
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        basis = A0_stack(:,:,i,j);
        norms(i,j) = sum(basis(:));
    end
end
norms = norms./sum(norms(:));
x_init = zeros(size(A0ft_stack));

% FISTA parameters
params.zeroMask = [];
params.zeroPad = [];
params.stoppingCriterion = 1;
params.tolerance = 1e-8;
params.L = 1;
params.lambda = 0.01;
params.beta = 1.1;
params.maxIter = 1000;
params.isNonnegative = 1;
params.noBacktrack = 0;
params.plotProgress = 0;

%% Create spots
% Produce 5x5 array of ring images
synth_spots = cell(P.num_var_t,P.num_var_r);
spot_fits = cell(P.num_var_t,P.num_var_r);
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        fprintf('Creating spots %i, %i\n',i,j);
        % Generate diffraction spot parameters
        P.theta_mean = round(P.num_theta/2);
        P.rad_mean = round(P.num_rad/2);

        % Generate image
        B = gaussian_basis_wrap_2D( P.num_theta,P.dtheta,P.theta_mean,P.var_theta(i),...
                                    P.num_rad,  P.drad,  P.rad_mean,  P.var_rad(j));
        % Add noise
        %B = B + randn(num_rad,num_theta)*mean(amplitudes)/10;

        synth_spots{i,j} = B;
        
        % Plot image
        % figure(2)
        % imshow(B,'DisplayRange',[0 300],'Colormap',jet)
    end
end

%% Fit each individual spot
spot_fits = cell(P.num_var_t,P.num_var_r);
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        fprintf('Fitting %i, %i\n',i,j);
        b = synth_spots{i,j};
        b = b/norm(b(:));
        [x_hat, err, obj, l_0]  = FISTA_Circulant(A0ft_stack,b,x_init,params);
        spot_fits{i,j} = x_hat;
        spots_l0(i,j) = l_0(end);
        spots_err(i,j) = err(end);
    end
end

%% Compute AWMV
true_awmv_az = zeros(P.num_var_t,P.num_var_r);
true_awmv_rad = zeros(P.num_var_t,P.num_var_r);
spot_awmv_az = zeros(P.num_var_t,P.num_var_r);
spot_awmv_rad = zeros(P.num_var_t,P.num_var_r);
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        fprintf('Computing AWMV %i, %i\n',i,j);
        xhat = spot_fits{i,j};
        [awmv_az, awmv_rad] = computeAWMV(xhat,P.var_theta,P.var_rad);
        spot_awmv_az(i,j) = awmv_az;
        spot_awmv_rad(i,j) = awmv_rad;
        true_awmv_az(i,j) = P.var_theta(i);
        true_awmv_rad(i,j) = P.var_rad(j);
%         fit_error
    end
end


%% Plot spots + fits 
figure(1)
for i = 1:P.num_var_t
    for j = 1:P.num_var_r
        spot = synth_spots{i,j};
        spot = spot/norm(spot(:));
        subplot(1,2,1)
        imshow(spot,'DisplayRange',[min(spot(:)) max(spot(:))],...
               'Colormap',jet,...
               'InitialMagnification','fit');
        xhat = spot_fits{i,j};
        spot_fit = Ax_ft_2D(A0ft_stack,xhat);
        subplot(1,2,2)
        imshow(spot_fit,'DisplayRange',[min(spot(:)) max(spot(:))],...
               'Colormap',jet,...
               'InitialMagnification','fit');
        vdf = squeeze(sum(sum(spot_fits{i,j},1),2))
        pause
    end
end