function D = distance_map(jTrue, nTrue, J, N)
% distance_map
%   Compute Euclidean distance from (jTrue,nTrue) to every pixel in JxN grid.
%
% Inputs:
%   jTrue : true scale index (1..J)
%   nTrue : true shift index (1..N)
%   J     : number of scales
%   N     : number of shifts
%
% Output:
%   D     : [J x N] distance matrix

    [jGrid, nGrid] = ndgrid(1:J, 1:N);   % grid of all coords
    D = sqrt((jGrid - jTrue).^2 + (nGrid - nTrue).^2);
end
