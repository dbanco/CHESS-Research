function [yn,y,N,M,T,Xtrue,Dtrue] = voigt_example_6peaks(sigma,plotFlag)
if nargin < 1
    sigma = 0.01;
end
if nargin < 2
    plotFlag = false;
end

rng(1234);

% General params
T = 50;
N = 105; M = 35;
center = (M+1)/2;
num_spots = 6;
K = num_spots; % For clarity

% --- Position and width trajectories ---
position = zeros(num_spots,T);
width = zeros(num_spots,T);
mixing = rand(num_spots,1);
posRange = linspace(15,90,num_spots);  % Base starting positions
y = zeros(N,T);

for k = 1:num_spots
    t_norm = linspace(0,1,T);
    % --- Position: slightly unique sinusoidal shift pattern ---
    shift_coef0 = 9*randn + posRange(k);
    shift_coef1 = 2 + randn;
    shift_coef2 = 8 + 3*randn;
    shift_pattern = shift_coef2*(t_norm).^2 + shift_coef1*t_norm + shift_coef0;
    position(k,:) = round(shift_pattern);

    % --- Width: polynomial (quadratic) growth ---
    width_coef0 = 2 + randn; % Add randomness
    width_coef2 = 4 + 2*randn; % Add randomness
    width_coef3 = 6 + 3*randn; % Add randomness
    w = width_coef0  + width_coef2*t_norm.^2 + width_coef3*t_norm.^3;
    
    width(k,:) = w;
    
    for t = 1:T
        y(:,t) = y(:,t) + voigt_basis_wrap_1D(N,position(k,t),width(k,t),mixing(k),'2-norm');
    end
end
yn = y + randn(N,T)*sigma;

% --- Plots ---
if plotFlag
    figure(1); imagesc(yn); title('Noisy signal');

    figure(2);
    num_slices = 8;
    times = round(linspace(1, T, num_slices));
    for i = 1:num_slices
        subplot(num_slices,1,i);
        plot(yn(:,times(i)), '-o');
        title(sprintf('Time slice %d', times(i)));
    end
end
Xtrue = [];
Dtrue = [];

end
