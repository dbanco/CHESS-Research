function visualizeTrajectories(X, K, J, threshold)
% X: [1 x N x KJ x T]
% K: number of atoms
% J: scales per atom
% threshold: energy threshold for blob detection

[~, N, KJ, T] = size(X);
colors = lines(K); % distinct color per atom

figure; hold on;
xlabel('Shift index (n)');
ylabel('Scale index (j)');
title('Independent trajectories per atom in scale-shift space');
set(gca, 'YDir','normal');

for k = 1:K
    subplot(K,1,k)
    scale_idx = ((k-1)*J + 1):(k*J);
    Xk = squeeze(X(1,:,scale_idx,:));  % [N x J x T]

    trajs = {}; % store blobs trajectories per atom

    % Extract blobs and CoMs independently at each time
    blobs_t = cell(T,1);
    for t = 1:T
        E = abs(Xk(:,:,t)).^2;
        BW = (E > threshold)';
        CC = bwconncomp(BW);

        CoMs = zeros(CC.NumObjects, 2); % [scale, shift]
        for b = 1:CC.NumObjects
            idx = CC.PixelIdxList{b};
            [j_idx, n_idx] = ind2sub([J,N], idx);
            E_blob = E(sub2ind([N,J], n_idx, j_idx));
            CoMs(b,1) = sum(j_idx .* E_blob) / sum(E_blob); % scale
            CoMs(b,2) = sum(n_idx .* E_blob) / sum(E_blob); % shift
        end
        blobs_t{t} = CoMs;
    end

    % Build trajectories per blob by linking only across consecutive frames for this atom
    trajs = {};
    for t = 1:(T-1)
        if t == 1
            for b = 1:size(blobs_t{1},1)
                trajs{b} = nan(T,2);
                trajs{b}(1,:) = blobs_t{1}(b,:);
            end
        end

        % Match blobs_t{t} to blobs_t{t+1} via nearest neighbor matching
        if isempty(blobs_t{t}) || isempty(blobs_t{t+1})
            continue;  % nothing to match, skip to next time step
        end
        distMat = pdist2(blobs_t{t}, blobs_t{t+1});
        for b = 1:length(trajs)
            [minDist, minIdx] = min(distMat(b,:));
            if ~isempty(minIdx) && minDist < Inf
                trajs{b}(t+1,:) = blobs_t{t+1}(minIdx,:);
                % Remove matched column so no double assignment
                distMat(:,minIdx) = Inf;
            else
                trajs{b}(t+1,:) = [NaN NaN];
            end
        end
    end

    % Plot each trajectory for this atom
    for b = 1:length(trajs)
        traj = trajs{b};
        valid = ~isnan(traj(:,1));
        plot(traj(valid,2), traj(valid,1), '-o', 'Color', colors(k,:), ...
            'DisplayName', sprintf('Atom %d Blob %d', k, b));
    end
    legend('Location','best');
    hold off;
end


end
