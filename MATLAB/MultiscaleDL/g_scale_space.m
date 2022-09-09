%% Construct gaussian scale space
close all
K = 5;
P.Ny = 16;
P.Nx = 64;
std1 = 0.3;
a1 = 1.5;
a2 = 1.8;
P.stds = [std1 std1;...
          std1*exp(1) std1*exp(1);...
          std1*exp(2) std1*exp(2);...
          std1*exp(3) std1*exp(3);...
          std1*exp(4) std1*exp(4)];
% P.stds = [std1 std1;...
%           std1*2^(1) std1*2^(1);...
%           std1*2^(2) std1*2^(2);...
%           std1*2^(3) std1*2^(3);...
%           std1*2^(4) std1*2^(4)];
% P.stds = [std1 std1;...
%           std1*a1^(1) std1*a2^(1);...
%           std1*a1^(2) std1*a2^(2);...
%           std1*a1^(3) std1*a2^(3);...
%           std1*a1^(4) std1*a2^(4)];
P.means = [P.Ny/2+1 P.Nx/2+1;...
           P.Ny/2+1 P.Nx/2+1;...
           P.Ny/2+1 P.Nx/2+1;...
           P.Ny/2+1 P.Nx/2+1;...
           P.Ny/2+1 P.Nx/2+1]; 
D = dictionary2D(P);

figure
for i = 1:K
    subplot(K,1,i)
    imagesc(D(:,:,i))
end

%% Convolve with gaussian kernel to get discrete scale space of each
P.Ns = 4;
Ds = zeros([size(D),P.Ns]);
Ds(:,:,1,1) = D(:,:,1);
Ds(:,:,2,1) = D(:,:,2);
figure

for i = 1:K
    for ss = 2:P.Ns
        Ds(:,:,i,ss) = conv2(Ds(:,:,i,ss-1),Ds(:,:,i,ss-1),'same');
    end
end

j = 1;
for i = 1:K
    for ss = 1:P.Ns
        subplot(K,P.Ns,j)
        imagesc(Ds(:,:,i,ss))
        j = j+1;
    end
end


