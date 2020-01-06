function x = forceMaskToZero( x, zeroMask )
%forceMaskToZero Forces zeroMask region to 0
if min(size(x)) > 1
    if(max(size(zeroMask))>0)
        x(zeroMask(:,1),zeroMask(:,2)) = 0;
    end
else
    if(max(size(zeroMask))>0)
        x(zeroMask) = 0;
    end 
end

end

