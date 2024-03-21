function showX(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% [N1,N2,KJ,T] = size(X);
% 
% K = 2;
% J = 8;
% 
% bigImg = 0.2*ones(N2*4,T*4);
% r1 = 1;
% r2 = N2;
% c1 = 1;
% c2 = T;
% for kj = 1:KJ  
%     bigImg(r1:r2,c1:c2) = X(1,:,kj,:);
%     c1 = c1 + T + 1;
%     c2 = c2 + T + 1;
%     if c2 > T*4+4
%         c1 = 1;
%         c2 = T;
%         r1 = r1 + N2 + 1;
%         r2 = r2 + N2 + 1;
%     end
% 
%    
% end
% 
% imagesc(bigImg)

vdf = squeeze(sum(squeeze(X),1));
imagesc(vdf)

%%
% xx = zeros(1,T);
% yy = xx;
% for t = 1:T
%     [x,y] = find(squeeze(X(:,:,1:8,t)));
%     xx(t) = x;
%     yy(t) = y;
% end
% 
% plot3(1:T,xx,yy,'-o')
% hold on
% 
% xx = zeros(1,T);
% yy = xx;
% for t = 1:T
%     [x,y] = find(squeeze(X(:,:,9:16,t)));
%     xx(t) = x;
%     yy(t) = y;
% end
% 
% plot3(1:T,xx,yy,'-o')
% 
% zlabel('\eta')
% ylabel('\sigma')
% xlabel('t')

%%
% figure
% hold on
% [r,c,v] = ind2sub([N2,8,T],find(squeeze(X(:,:,1:8,:)>1e-4)) );
% plot3(v,r,c,'o')
% 
% [r,c,v] = ind2sub([N2,8,T],find(squeeze(X(:,:,9:16,:)>1e-4)) );
% plot3(v,r,c,'o')
% ylabel('\eta')
% zlabel('\sigma')
% xlabel('t')


end