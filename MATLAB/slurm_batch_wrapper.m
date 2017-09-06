function slurm_batch_wrapper( index )
%batch_wrapper Runs function as specified by generic inout files

index_str = num2str(index);

% Load inputs
load(['varin_',index_str,'.mat'])

% Evaluate function
funcH = eval(['@',funcName]);
funcH(varin{1},varin{2},varin{3},varin{4})

end

