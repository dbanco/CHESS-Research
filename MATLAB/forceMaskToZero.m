function x = forceMaskToZero( x, zeroMask )
%forceMaskToZero Forces zeroMask region to 0
try
    if(max(size(zeroMask))>0)
        x(zeroMask(:,1),zeroMask(:,2)) = 0;
    end
catch
    if(max(size(zeroMask))>0)
        x(zeroMask(:,2)) = 0;
    end 
end

end

