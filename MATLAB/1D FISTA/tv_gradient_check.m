close all
N = 100;
t = 20;
deltaX = 1/t;
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
figure(1)
plot(f')
D = [-1 1 0;...
      0 -1 1]/deltaX;
total = sum(x2(:));

figure(2)
plot( sum(((D*f)').^2,2) )


gradJ = gradientTV(f,deltaX,tvBeta,D);
gradTV = zeros(size(gradJ));
for j = 1:t
    for k = 1:t
        if j == k
            gradTV(j) = gradTV(j) + gradJ(k)*( 1/total - vdf2(j)/total );
        else
            gradTV(j) = gradTV(j) - gradJ(k)*vdf2(j)/total;
        end     
    end
end



%% finite difference check 1

gradF = zeros(size(gradJ));
for j = 1:t

        vdf22 = vdf2;
        vdf22(j) = vdf22(j) + 1e-8;
        f22 = [vdf1;vdf22;vdf3];
        J = deltaX*sum(sqrt( (D*f).^2 + tvBeta^2));      
        J22 = deltaX*sum(sqrt( (D*f22).^2 + tvBeta^2));
        gradF(j) = (J(j)-J22(j))/1e-8;
end

gradJ-gradF

%% finite difference check 2
% close all
epsi = 1e-8;
gradF2 = zeros(size(gradJ'));

for i = 1:N
	for j = 1:t
        x22 = x2;
        x22(i,j) = x22(i,j) + epsi;
        vdf22 = sum(x22,1)./sum(x22(:));
        f22 = [vdf1;vdf22;vdf3];
        J = deltaX*sum(sqrt( (D*f).^2 + tvBeta^2));      
        J22 = deltaX*sum(sqrt( (D*f22).^2 + tvBeta^2));
        gradF2(j) = gradF2(j) + (J(j)-J22(j))/epsi;
    end
end
N*gradTV-gradF2'


figure(3)

subplot(2,1,1)
hold on
plot(gradJ,'Linewidth',4)
plot(gradF,'Linewidth',2)
subplot(2,1,2)
hold on
plot(N*gradTV,'Linewidth',4)
plot(gradF2','Linewidth',2)