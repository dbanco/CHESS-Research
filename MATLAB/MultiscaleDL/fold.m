function [X] = fold( T, dim, i )
X = shiftdim(reshape(T,circshift(dim,1-i)),numel(dim)-i+1);
end