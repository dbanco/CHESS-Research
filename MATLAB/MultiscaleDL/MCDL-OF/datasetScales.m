function [K,J,scales] = datasetScales(dataset)

switch dataset
    case 'dissertation_adjust2'
        K = 2;
        scales = cell(K,1);
        scales{1} = genRationals([0;1],[1;1],8,8, 1/6);
        scales{2} = genRationals([0;1],[1;1],8,8, 1/6);
        J = size(scales{1},2);
    case 'pseudo-voigt_unmatched'
        K = 1;
        scales = cell(K,1);
        scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
        J = size(scales{1},2);
end

end