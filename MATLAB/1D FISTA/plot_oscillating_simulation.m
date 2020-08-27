


filedir = 'E:\CHESS_data\simulated_two_spot_1D_gnoise6_osc_1';
base_file = 'polar_vector_%i.mat';

osc_fig = figure(1);
[ha1, pos1] = tight_subplot(5,20,[0.1 0.03],[.02 .08],[.02 .02]); 

for i = 1:100
    fname = sprintf(base_file,i);
    load(fullfile(filedir,fname))
    axes(ha1(i))
    plot(polar_vector)
end