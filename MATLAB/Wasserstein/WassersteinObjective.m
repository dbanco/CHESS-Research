function [wObj, T] = WassersteinObjective( vdf, vdfs, lam, D )
%WassersteinObjective Summary of this function goes here
%   Detailed explanation goes here
wObj = 0;

for i = 1:numel(vdfs)
    neighbor_vdf = vdfs{i};
    [ ~, Wd, T ] = WassersteinGrad(vdf, neighbor_vdf(:), lam, D);
    wObj = wObj + Wd; 
end

end

