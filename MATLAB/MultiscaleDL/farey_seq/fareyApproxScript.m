%% Farey approximation
% Input number 1-3
% Output nearest rational number with numerator ~<100
N = 10000;
testM = 10:200;
M = numel(testM);
error = zeros(N,M);
for m = 1:M
    for n = 1:N
        input = 2*rand(1)+1;
        [c1,c2] = fareyApprox(input,m);
        error(n,m) = abs(input-c1/c2);
    end
end

loglog(testM,mean(error))
hold on
loglog(testM,std(error))
legend('avg error','std error')
xlabel('denominator limit')
title('10^4 trials of rational approximations of x \in [1,3]')
grid on