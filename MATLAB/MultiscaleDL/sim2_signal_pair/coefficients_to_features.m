function features = coefficients_to_features(outputs)
% Convert sparse coefficients in output struct to feature vectors
% Each nonzero coefficient becomes a feature vector:
% [amplitude, shift_index, scale_index, time_index, atom_index]

N  = outputs.N;
M  = outputs.M;
K  = outputs.K;
J  = outputs.J;
X  = squeeze(outputs.X); % (N+M-1) x (K*J) x T

[~, KJ, T] = size(X);

% Preallocate (worst case all nonzero)
nzmax = nnz(X);
features = zeros(nzmax,5,T);

count = 0;
for t = 1:T
    for kj = 1:KJ
        [shift_idx, ~, amp] = find(X(:,kj,t));
        if ~isempty(shift_idx)
            for ii = 1:numel(shift_idx)
                count = count + 1;
                % Map kj back to (atom k, scale j)
                atom_idx  = ceil(kj / J);
                scale_idx = mod(kj-1, J) + 1;
                features(count,:,t) = [ ...
                    amp(ii), ...        % amplitude
                    shift_idx(ii), ...  % shift index
                    scale_idx, ...      % scale index
                    t, ...              % time index
                    atom_idx];          % atom index
            end
        end
    end
end

% Trim unused rows
features = features(1:count,:,:);

end