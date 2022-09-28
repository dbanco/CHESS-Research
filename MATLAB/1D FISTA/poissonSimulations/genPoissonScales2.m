function [alpha_val,theta_stds] = genPoissonScales2(N,T,levels,sim)
%% Comptue error for different data scales
Reps = 5;
alpha_val = zeros(numel(levels),1);
for  i= 1:numel(levels)
    a = 0.1;
    b = 10000;
    snrAvgA = 0;
    snrAvgB = 0;
    snrAvgC = 0;
    for rep = 1:Reps
        [~,~,theta_stds,~,~,snr] = genSimDataPoisson2(N,T,a,sim);
        snrAvgA = snrAvgA + mean(snr)/Reps;
        [~,~,theta_stds,~,~,snr] = genSimDataPoisson2(N,T,b,sim);
        snrAvgB = snrAvgB + mean(snr)/Reps;
        c = (a+b)/2;
        [~,~,theta_stds,~,~,snr] = genSimDataPoisson2(N,T,c,sim);
        snrAvgC = snrAvgC + mean(snr)/Reps;
    end
    crit = abs(snrAvgC-levels(i));
        
    while crit>1e-2
        if sign(snrAvgA-levels(i)) == sign(snrAvgC-levels(i))
            a = c;
            snrAvgA = snrAvgC;
        elseif sign(snrAvgB-levels(i)) == sign(snrAvgC-levels(i))
            b = c;  
            snrAvgB = snrAvgC;
        else
            error('Bisection failed')
        end
        c = (a+b)/2;
        snrAvgC = 0;
        for rep = 1:Reps
            [~,~,theta_stds,~,~,snr] = genSimDataPoisson2(N,T,c,sim);
            snrAvgC = snrAvgC + mean(snr)/Reps;  
        end      
        crit = abs(snrAvgC-levels(i));
    end
    alpha_val(i) = c;
end