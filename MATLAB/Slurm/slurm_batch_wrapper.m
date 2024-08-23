function slurm_batch_wrapper( index )
%batch_wrapper Runs function as specified by generic inout files

index_str = num2str(index);

% Load inputs
load(['varin_',index_str,'.mat'])

% Evaluate function
funcH = eval(['@',funcName]);
switch numel(varin)
    case 1
        funcH(varin{1})
    case 2
        funcH(varin{1},varin{2})
    case 3
        funcH(varin{1},varin{2},varin{3})
    case 12
        funcH(varin{1},varin{2},varin{3},varin{4},varin{5},varin{6},...
            varin{7},varin{8},varin{9},varin{10},varin{11},varin{12})
    case 16
        funcH(varin{1},varin{2},varin{3},varin{4},varin{5},varin{6},...
            varin{7},varin{8},varin{9},varin{10},varin{11},varin{12},...
            varin{13},varin{14},varin{15},varin{16})
end