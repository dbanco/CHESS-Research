function full_array = fill_array( full_array, idx_array, idx )
%fill_array Updates full array at indices specified by idx from idx_array
for k = 1:numel(idx)
    [i,j] = ind2sub(size(full_array),idx(k));
    full_array(i,j) = idx_array(k);
end

