function D = grid_distance_matrix(J,N)
% grid_distance_matrix
%   Compute full (JN x JN) distance matrix
%   between all pixels in a JxN grid.

    % Coordinates of each pixel (flattened)
    [jGrid,nGrid] = ndgrid(1:J,1:N);
    coords = [jGrid(:), nGrid(:)];   % [JN x 2]

    % Compute pairwise distances
    D = pdist2(coords, coords);      % [JN x JN]
end
