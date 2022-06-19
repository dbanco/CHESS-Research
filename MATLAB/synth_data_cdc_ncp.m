%% Spot example
% Ring sampling parameters
N1 = 16*4;
N2 = 100*4;
T = 50;

% Basis function variance parameters
num_spots = 3;
theta1 = linspace(1,0,T)';
theta = [theta1.^2 (1-theta1.^2)];
std1   = [0.52 0.56 0.5;
          0.82 0.61 0.51]*4;
std2 = [2, 3, 1.5;
        8, 6, 9;]*4;
mean1 = [3,4,2;
        11,12,13]*4;
mean2 = [15, 45, 65;
         35, 65, 85]*4;

std1_t = theta*std1;
std2_t = theta*std2;
mean1_t = theta*mean1;
mean2_t = theta*mean2;

B = zeros(N1,N2,T);
for t = 1:T
    for i = 1:num_spots
        B(:,:,t) = B(:,:,t) +...
        gaussian_basis_wrap_2D( N1,mean1_t(t,i),std1_t(t,i),...
                                N2,mean2_t(t,i),std2_t(t,i),'rms');  
    end
%     B(:,:,t)= B(:,:,t)+B(:,:,t).*randn(N1,N2)/20; 
end

f = figure(2);
for t = 1:T
    imagesc(B(:,:,t))
    f.Position = [100 100 540 300];
    pause(0.05)
end