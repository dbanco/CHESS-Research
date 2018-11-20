function wObj = WassersteinObjective( vdf, vdfs, lam, D )
%WassersteinObjective Summary of this function goes here
%   Detailed explanation goes here
wObj = 0;

for i = 1:numel(vdfs)
    neighbor_vdf = vdfs{i};
    Wd = sinkhornKnoppTransport(D, lam, vdf, neighbor_vdf(:));
    wObj = wObj + Wd; 
end

end

