function [overlapMinRow,overlapMaxRow,overlapMinCol,overlapMaxCol] = overlapMatrix(A0ft_stack,T,R)
%overlapMatrix Summary of this function goes here

% Output matrices
overlapMinRow = zeros(T,R,T,R);
overlapMaxRow = zeros(T,R,T,R);
overlapMinCol = zeros(T,R,T,R);
overlapMaxCol = zeros(T,R,T,R);

for t = 1:T
    for r = 1:R
        for t2 = 1:T
            for r2 = 1:R
            % Convolve atoms to identify overlap
            atom_ft_1 = squeeze(A0ft_stack(:,:,t,r));
            atom_ft_2 = squeeze(A0ft_stack(:,:,t2,r2));
            conv_atom = real(ifft2(atom_ft_1.*atom_ft_2));
            [Drow,Dcol] = find(conv_atom > 1e-8); 
            overlapMinRow(t,r,t2,r2) = min(Drow);
            overlapMaxRow(t,r,t2,r2) = max(Drow);
            overlapMinCol(t,r,t2,r2) = min(Dcol);
            overlapMaxCol(t,r,t2,r2) = max(Dcol);
            end
        end
    end
end

end

