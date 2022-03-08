%% Load MMPAD sequence
r = 1;
data_dir = ['D:\MMPAD_data_nr1\ring', num2str(r), '_zero'];
T = 546;

for t = 1:T
	load(fullfile(data_dir,['mmpad_img_',num2str(t),'.mat']))
    if t == 1
        [N,M] = size(polar_image);
        X = zeros(N,M,T);
    end
    X(:,:,t) = polar_image;
end

%% Tucker ALS compressibility
R1 = 12; 
R2 = 60; 
R3 = 100;

R = 100;
X = tensor(X);
% rel_error = zeros(R,1);
% ranks = zeros(R,3);
% t_rank = zeros(R,1);
for i = 1:300
    if R1 == N
        if R2 == M
            R3 = R3 + 1;
        else
            switch mod(i,2)
                case 1
                    R2 = R2 + 1;
                otherwise
                    R3 = R3 + 1;
            end
        end
    else 
        switch mod(i,8)
            case 0 
                R1 = R1 + 1;
            case 2
                R2 = R2 + 1;
            case 5
                R2 = R2 + 1;
            otherwise
                R3 = R3 + 1;
        end
    end
            Tuck = tucker_als(X,[R1,R2,R3]);
            rel_error(i) = norm(X-Tuck)/norm(X);
            ranks(i,:) = [R1, R2, R3];
            t_rank(i) = R1*R2*R3;
end

plot(t_rank,rel_error,'o-')
xlabel('Rank')
ylabel('Error')

%% Load MMPAD sequence
r = 1;
data_dir = ['D:\MMPAD_data_nr1\ring', num2str(r), '_zero'];
T = 546;
X_norm = X;
X_sum = X;
for t = 1:T
    X_norm(:,:,t) = X(:,:,t)/sqrt(sum(X(:,:,t).^2,'all'));
    X_sum(:,:,t) = X(:,:,t)/sum(X(:,:,t),'all');
end
X_all = zeros(N,M,T,3);
X_all(:,:,:,1) = X;
X_all(:,:,:,2) = X_norm;
X_all(:,:,:,3) = X_sum;

%% Compute Tucker decomposition
R1 = 10; R2 = 20; R3 = 30;
time_err = zeros(T,3);
for k = 1:3
    X_t = tensor(squeeze(X_all(:,:,:,k)));
    % Tucker Decomposition and compare
    
    Tuck = tucker_als(X_t,[R1,R2,R3]);
    x = X_t.data;
    tuck = tensor(Tuck);
    tuck = tuck.data;
    % Plot a frame
   
    for i = 1:T
        time_err(i,k) = norm(x(:,:,i)-tuck(:,:,i))/norm(x(:,:,i));
    end
%     figure(1)
%     for i = 60
%         frame1 = [x(:,:,i);tuck(:,:,i)];
%         imagesc(frame1)
%         title(sprintf('Frame %i',i))
%     end
end
figure(2)
plot(time_err(:,3))
xlabel('time')
ylabel('Reconstruction error')
legend(sprintf('natural: %0.2f',mean(time_err(:,1))),...
        sprintf('2-norm: %0.2f',mean(time_err(:,2))),....
        sprintf('1-norm: %0.2f',mean(time_err(:,3))))
    
%% Inspect Tucker Decomp 
R1 = 10; R2 = 20; R3 = 20;
k = 2;
X_t = tensor(squeeze(X_all(:,:,:,k)));
% Tucker Decomposition and compare

Tuck = tucker_als(X_t,[R1,R2,R3]);

% Visualize factors
figure(3)
for j = 1:3
    subplot(1,3,j)
    imagesc(Tuck.U{j})
end

% Compute error
x = X_t.data;
tuck = tensor(Tuck);
tuck = tuck.data;
for i = 1:T
    time_err(i,k) = norm(x(:,:,i)-tuck(:,:,i))/norm(x(:,:,i));
end
mean(time_err(:,2))

% View a frame
    figure(1)
    for i = 60
        frame1 = [x(:,:,i);tuck(:,:,i)];
        imagesc(frame1)
        title(sprintf('Frame %i',i))
    end
%% T-SVD compressibility

