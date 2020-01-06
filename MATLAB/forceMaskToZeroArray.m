function x = forceMaskToZeroArray( x, zeroMask )
%forceMaskToZero Forces zeroMask region to 0

if numel(size(x)) == 4
    if(max(size(zeroMask))>0)
        x(zeroMask(:,1),zeroMask(:,2),:,:) = 0;
    end
else
    if(max(size(zeroMask))>0)
        x(zeroMask,:) = 0;
    end
end


end

