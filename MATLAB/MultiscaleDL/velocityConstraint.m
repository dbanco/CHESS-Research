% Compute constraint to follow velocities
% True solution
load('C:\Users\dpqb1\Documents\Outputs\toy1_exp_OF1vel1_matched_200exact_sig_1\output_j1_sig_0.00e+00_lam1_0.00e+00_lam2_0.00e+00.mat');
data = squeeze(outputs.X);
[N, KJ, T] = size(data);
[u,v,Fy,Fx,Ft]  = computeHornSchunkDict(data,K,smoothness,maxIt);
K = 2;

velCon(data,K)

%% Show
t = 3;
figure;
imagesc(Vt1(:,:,t))

figure;
imagesc(Vtm1(:,:,t))

figure;
imagesc(data(:,:,t))

figure;
imagesc(data(:,:,t+1))


figure;
imagesc( Vt1(:,:,t).*data(:,:,t+1) + Vtm1(:,:,t).*data(:,:,t) )

%% Place -1 where velocities are defined for Vtm1
% Vtm1(abs(u(:,:,2:T))>0) = -1;
% Vtm1(abs(v(:,:,2:T))>0) = -1;
%     
% % Place 1 where velocities point towards
% for i = 1:N
%     for j = 1:KJ
%         for t = 2:T
%             % Check x-direction
%             if u(i,j,t) > 0
%                 jj = j+1;
%             elseif u(i,j,t) < 0
%                 jj = j-1;
%             else
%                 jj = j;
%             end
%             % Check y-direction
%             if v(i,j,t) > 0
%                 ii = i+1;
%             elseif v(i,j,t) < 0
%                 ii = i-1;
%             else
%                 ii = i;
%             end
%             % Set value
%             Vt1(ii,jj,t-1) = Vt1(ii,jj,t-1) + 1;
%         end
%     end
% end
% 
% % Define where velocitys are not to penalize more heavily
% for t = 2:T
%     [is,js] = find((Vtm1(:,:,t-1)>0));
%     noV(is,js,t-1) = 1;
%     [is,js] = find((Vt1(:,:,t-1)>0));
%     noV(is,js,t) = 1;
% end
% 
% % Compute velocity constraint
% sum((Vt1.*data(:,:,2:T) - Vtm1.*data(:,:,1:T-1)).^2,'all')

