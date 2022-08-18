%% Show simulated data
filedir = 'E:\CHESS_data\simulated_two_spot_1D_gnoise6_osc_1';
base_file = 'polar_vector_%i.mat';

i = 1;
fname = sprintf(base_file,i);
load(fullfile(filedir,fname))

% 1st image
figure(1)
plot(polar_vector)
ylabel('Intensity')
xlabel('\eta')

filedir = 'E:\CHESS_data\simulated_two_spot_1D_gnoise6_osc_11';
true_data = load(fullfile(filedir,'synth_data.mat'))
true_awmv = zeros(150,1);
for j = 1:150
    true_stds = true_data.synth_sample{j}.std_theta;
    true_amps = true_data.synth_sample{j}.amplitudes;
    true_awmv(j) = true_stds'*true_amps/sum(true_amps);
end

% True AWMV
figure(2)
plot(true_awmv)

