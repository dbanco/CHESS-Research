function img_pad = zeroPad( img, zeroPad )
%zeroPad Zero pads image

if size(zeroPad) > 1
    [n,m] = size(img);
    if(max(zeroPad > 0))
        % Pad image
        img_pad = zeros(2*zeroPad(1)+n,2*zeroPad(2)+m);
        img_pad(zeroPad(1)+1:zeroPad(1)+n,zeroPad(2)+1:zeroPad(2)+m) = img;
    else
        % Do not pad
        img_pad = img;
    end
else
    n = numel(img);
    if(max(zeroPad > 0))
        % Pad image
        img_pad = zeros(2*zeroPad+n);
        img_pad(zeroPad+1:zeroPad+n) = img;
    else
        % Do not pad
        img_pad = img;
    end  
end

end

