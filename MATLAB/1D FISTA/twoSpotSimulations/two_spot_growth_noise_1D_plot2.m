%% Load spread data
clear all

num_imgs = 10;
spreadDir = fullfile('D:','CHESS_data','spread_results');

datasetArray = {'simulated_two_spot_1D_noise_independent',...
                'simulated_two_spot_1D_noise_independent_quadratic',...
                'noise_1D_3_1a',...
                'noise_1D_3_2a'};

line = [3,6,1,1];

noise_std = [0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3];

legend_str = {'Unreg','Unreg Quad','Reg 0.25','Reg 0.5','Truth'};

n = numel(datasetArray);
colors = jet(numel(legend_str));

for i = 1:n
    spread_dir = fullfile(spreadDir,['spread_',datasetArray{i},'.mat'])
    ring_data{i} = load(spread_dir);
end

%% Load truth
load(['D:\CHESS_data\simulated_two_spot_1D_noise_12\synth_data.mat'])
truth_awmv_az = zeros(num_imgs,1);

for i = 1:num_imgs
    sample = synth_sample{i};
    az_awmv = sample.std_theta'*sample.amplitudes/sum(sample.amplitudes(:));
    truth_awmv_az(i) = az_awmv;
end

%% Plot time varying spread and error functions
cv_i = 1;
for i = 1:n
    nn = numel(noise_std);
    for ii = line(i)     
        figure(2)
        cv = colors(cv_i,:);
        plot(ring_data{i}.az_spread(:,ii),'-o','Color',cv,'MarkerSize',3)
        hold on
        title('Azimuthal AWMV')
        xlabel('Time')

        figure(3)
        cv = colors(cv_i,:);
        plot(ring_data{i}.rel_err(:,ii),'-o','Color',cv,'MarkerSize',3)
        hold on
        title('Relative Error')
        xlabel('Time')

        figure(4)
        cv = colors(cv_i,:);
        semilogy(ring_data{i}.sparsity(:,ii),'-o','Color',cv)
        hold on
        title('Sparsity')
        xlabel('Time')
        
        figure(5)
        cv = colors(cv_i,:);
        semilogy(ring_data{i}.objective(:,ii),'-o','Color',cv)
        hold on
        title('Objective')
        xlabel('Time')
        
        cv_i = cv_i + 1;
    end
end
%% Truth lines
    
figure(2)
plot(truth_awmv_az,'-x','MarkerSize',2)
legend(legend_str,'Location','Best')

figure(3)
legend(legend_str,'Location','Best')

figure(4)
legend(legend_str,'Location','Best')

figure(5)
legend(legend_str,'Location','Best')
% figure(6)
% for j = 1
%     plot(ring_data{j}.wass_dist,'-')
%     hold on
% end
% legend('1','2','3','4','Location','Best')
% title('Wasserstein distance to neighbor')
% xlabel('Time')