close all
N = 100;
t = 20;
tvBeta = 1e-4;
t0 = 2;
a = 1:t;
x1 = repmat(a,N,1);
x2 = repmat((a-(t+1)/2).^2,N,1);
x3 = repmat(t+1-a,N,1);
vdf1 = sum(x1,1)./sum(x1(:));
vdf2 = sum(x2,1)./sum(x2(:));
vdf3 = sum(x3,1)./sum(x3(:));

f = [vdf1;vdf2;vdf3];
fx = [x1(:), x2(:), x3(:)]';
figure(1)
plot(f')
D = [-1 1 0;...
      0 -1 1];
total = sum(x2(:));

figure(2)
plot( sum(((D*f)').^2,2) )


% gradJ = gradientTV(f,deltaX,tvBeta,D);
% gradTV = zeros(size(gradJ));
% for j = 1:t
%     for k = 1:t
%         if j == k
%             gradTV(j) = gradTV(j) + gradJ(k)*( 1/total - vdf2(j)/total );
%         else
%             gradTV(j) = gradTV(j) - gradJ(k)*vdf2(j)/total;
%         end     
%     end
% end
% 


%% finite difference check 1

% gradF = zeros(size(gradJ));
% for j = 1:t
% 
%         vdf22 = vdf2;
%         vdf22(j) = vdf22(j) + 1e-8;
%         f22 = [vdf1;vdf22;vdf3];
%         J = deltaX*sum(sqrt( (D*f).^2 + tvBeta^2));      
%         J22 = deltaX*sum(sqrt( (D*f22).^2 + tvBeta^2));
%         gradF(j) = -(J(j)-J22(j))/1e-8;
% end
% 
% gradJ-gradF

%% finite difference check {vdf1i - vdf2i|
close all
epsi = 1e-8;
gradF2 = zeros(size(x1));

for i = 1:N
	for j = 1:t
        x22 = x2;
        x22(i,j) = x22(i,j) + epsi;
        vdf22 = sum(x22,1)./sum(x22(:));
        f22 = [vdf1;vdf22;vdf3];
        J = sum(sum(sqrt( (D*f).^2 + tvBeta^2)));      
        J22 = sum(sum(sqrt( (D*f22).^2 + tvBeta^2)));
        gradF2(i,j) = gradF2(i,j) - (J-J22)/epsi;
    end
end

gradTV = zeros(size(x1));
gradJ = gradientTV(f,tvBeta,D);
for i = 1:N
    for j = 1:t
        for k = 1:t
            if j == k
                gradTV(i,j) = gradTV(i,j) + gradJ(k)*( 1/total - vdf2(k)/total );
            else
                gradTV(i,j) = gradTV(i,j) - gradJ(k)*vdf2(k)/total;
            end 
        end
    end
end
abs(mean(gradTV(:) - gradF2(:)))

%% finite difference check |x1i-x2i|
% close all
epsi = 1e-8;
gradF3 = zeros(size(x1));

for i = 1:N
	for j = 1:t
        x22 = x2;
        x22(i,j) = x22(i,j) + epsi;
        fx22 = [x1(:), x22(:), x3(:)]';
        J = sum(sum(sqrt( (D*fx).^2 + tvBeta^2)));      
        J22 = sum(sum(sqrt( (D*fx22).^2 + tvBeta^2)));
        gradF3(i,j) = gradF3(i,j) - (J-J22)/epsi;
    end
end
gradTVx = gradientTV(fx,tvBeta,D);
gradTVx = reshape(gradTVx,size(x1));
abs(mean(gradTVx(:) - gradF3(:)))

figure(3)
subplot(2,1,1)
hold on
plot(gradTV(:),'Linewidth',4)
plot(gradF2(:),'Linewidth',2)
legend('TV','FD')

subplot(2,1,2)
hold on
plot(gradTVx(:),'Linewidth',4)
plot(gradF3(:),'Linewidth',2)
legend('TVx','FD')