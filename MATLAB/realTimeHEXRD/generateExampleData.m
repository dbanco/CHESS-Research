function [B,B_noise,awmv_true] = generateExampleData(N,T)
%generateExampleData Generate example poisson measurements of gaussian
%intensity peaks

numSpots = 2;
B = zeros(N,T);
B_noise = zeros(N,T);
amplitude = 80*[0.4 0.7]+1;
mean_param = N*[0.3 0.7];
widths = [5 8];

awmv_true = zeros(T,1);
for t = 1:T
    for i = 1:numSpots
        b = gaussian_basis_wrap_2D(N,...
                                   mean_param(i),...
                                   widths(i),...
                                   '2-norm');
       awmv_true(t) = awmv_true(t) + amplitude(i)*widths(i);
       B(:,t) = B(:,t) + amplitude(i)*b;
    end
    awmv_true(t) = awmv_true(t)/sum(amplitude(:));
    B_noise(:,t) = poissrnd(B(:,t));
end
end

