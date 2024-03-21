function [B,B_noise,awmv_true1,awmv_true2] = generateExampleData2D(N1,N2,T,numSpots)
%generateExampleData Generate example poisson measurements of gaussian
%intensity peaks

if nargin<4
    numSpots = 2;
end
B = zeros(N1,N2,T);
B_noise = zeros(N1,N2,T);
amplitude = 80*[0.4 0.7]+1;

mean_param1 = N1*[0.3 0.7];
widths1 = [5 8];
mean_param2 = N2*[0.3 0.7];
widths2 = [5 8];

awmv_true1 = zeros(T,1);
awmv_true2 = zeros(T,1);
for t = 1:T
    for i = 1:numSpots
        b = gaussian_basis_wrap_2D(N1,mean_param1(i),widths1(i),...
                                   N2,mean_param2(i),widths2(i),'2-norm');
       awmv_true1(t) = awmv_true1(t) + amplitude(i)*widths1(i);
       awmv_true2(t) = awmv_true2(t) + amplitude(i)*widths2(i);
       B(:,:,t) = B(:,:,t) + amplitude(i)*b;
    end
    awmv_true1(t) = awmv_true1(t)/sum(amplitude(:));
    awmv_true2(t) = awmv_true2(t)/sum(amplitude(:));
    B_noise(:,:,t) = poissrnd(B(:,:,t));
end
end

