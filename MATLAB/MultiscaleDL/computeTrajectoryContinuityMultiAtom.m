function objTraj = computeTrajectoryContinuityMultiAtom(X, u, v, K, J, threshold)
% X: [1 x N x KJ x T]
% u, v: same size
% K: number of atoms
% J: number of scale bins per atom
% threshold: scalar for energy threshold

[~, N, KJ, T] = size(X);

% Sanity check
if K*J ~= KJ
    error('K * J must equal KJ!');
end

objTraj = 0;

for k = 1:K
    % Indices for this atom's scale slices
    scale_idx = ((k-1)*J + 1):(k*J);

    % Extract [N x J x T]
    Xk = squeeze(X(1,:,scale_idx,:));  % [N x J x T]
    uk = squeeze(u(:,scale_idx,:));
    vk = squeeze(v(:,scale_idx,:));

    for t = 1:(T-1)
        % Energy maps: [N x J]
        E_now = abs(Xk(:,:,t)).^2;
        E_next = abs(Xk(:,:,t+1)).^2;

        % Permute to [J x N] for rows = scale, cols = shift
        BW_now = E_now' > threshold;  % logical [J x N]
        BW_next = E_next' > threshold;

        CC_now = bwconncomp(BW_now);
        CC_next = bwconncomp(BW_next);

        for b = 1:CC_now.NumObjects
            idx = CC_now.PixelIdxList{b};
            [j_idx, n_idx] = ind2sub([J,N], idx);

            % Extract blob energy: map back to [N x J] indexing
            E_blob = E_now(sub2ind([N,J], n_idx, j_idx));

            % CoM for this blob
            shiftCoM_now = sum(n_idx .* E_blob) / sum(E_blob);
            scaleCoM_now = sum(j_idx .* E_blob) / sum(E_blob);

            % Average flow
            u_blob = uk(:,:,t)';
            v_blob = vk(:,:,t)';

            avg_u = sum(u_blob(idx) .* E_blob) / sum(E_blob);
            avg_v = sum(v_blob(idx) .* E_blob) / sum(E_blob);

            % ---- Next frame ----
            best_dist = inf;
            shiftCoM_next = NaN;
            scaleCoM_next = NaN;

            for b2 = 1:CC_next.NumObjects
                idx2 = CC_next.PixelIdxList{b2};
                [j2_idx, n2_idx] = ind2sub([J,N], idx2);

                E_blob2 = E_next(sub2ind([N,J], n2_idx, j2_idx));

                shiftCoM2 = sum(n2_idx .* E_blob2) / sum(E_blob2);
                scaleCoM2 = sum(j2_idx .* E_blob2) / sum(E_blob2);

                dist = (shiftCoM2 - scaleCoM_now)^2 + (shiftCoM2 - shiftCoM_now)^2;

                if dist < best_dist
                    best_dist = dist;
                    shiftCoM_next = shiftCoM2;
                    scaleCoM_next = scaleCoM2;
                end
            end

            if ~isnan(shiftCoM_next)
                % Add penalty
                objTraj = objTraj + ...
                    (shiftCoM_next - shiftCoM_now - avg_v)^2 + ...
                    (scaleCoM_next - scaleCoM_now - avg_u)^2;
            end
        end
    end
end

end
