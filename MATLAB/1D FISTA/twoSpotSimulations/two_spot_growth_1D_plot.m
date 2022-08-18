%% Load spread data
clear all

num_imgs = 20;
spreadDir = fullfile('D:','CHESS_data','spread_results');

datasetArray = {'two_spot_growth_1D_25_renorm_1-6'...
                'two_spot_growth_1D_25_l2norm_1-4'};
lines = [6,4];
gamma_vals = {[0.15 0.1 0.08 0.05 0.02 0.01],...
              [0.15 0.1 0.08 0.05]};

legend_str = cell(sum(lines)+1,1);
k = 1;
for i = 1:numel(gamma_vals)
    g_vals = gamma_vals{i};
    for ii = 1:numel(g_vals)
        legend_str{k} = num2str(g_vals(ii));
        k = k + 1;
    end
end
legend_str{end} = 'truth';

n = numel(datasetArray);
nn = 5;
colors = jet(numel(legend_str));

for i = 1:n
    spread_dir = fullfile(spreadDir,['spread_',datasetArray{i},'.mat'])
    ring_data{i} = load(spread_dir);
end

%% Load truth
load(['D:\CHESS_data\simulated_data_two_spot_growth_1D_25\synth_data.mat'])
truth_awmv_az = zeros(num_imgs,1);

for i = 1:num_imgs
    sample = synth_sample{i};
    az_awmv = sample.std_theta'*sample.amplitudes/sum(sample.amplitudes(:));
    truth_awmv_az(i) = az_awmv;
end

%% Plot time varying spread and error functions
cv_i = 1;
for i = 1:n
    nn = lines(i);
    for ii = 1:nn     
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