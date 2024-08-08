function x = forcePadToZero( x, zeroPad )
%forcePadToZero Forces zeroPad region to 0

if(max(zeroPad)>0)
    x(1:zeroPad(1),:) = 0;
    x((end-zeroPad(1)+1):end,:) = 0;
    x(:,1:zeroPad(2)) = 0;
    x(:,(end-zeroPad(2)+1):end) = 0;
end

end

