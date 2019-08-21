%% Load spread data
clear all

num_imgs = 10;
spreadDir = fullfile('D:','CHESS_data','spread_results');

datasetArray = {'mmpad_1D_subset_1-4',...
                'mmpad_1D_subset2_1-4',...
                'mmpad_1D_subset_independent'}
lines = [4,4,1];
gamma_vals = {[0.1 0.08 0.05 0.01],...
              [10 5 1 0.5],...
              [0]};

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
%% Truth lines + legends
    
figure(2)
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