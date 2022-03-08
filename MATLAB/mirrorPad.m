function img_pad = mirrorPad( img, mPad )
%zeroPad Zero pads image

if numel(mPad) > 1
    [n,m] = size(img);
    if(max(mPad > 0))
        % Pad image
        img_pad = zeros(2*mPad(1)+n,2*mPad(2)+m);
        img_pad(mPad(1)+1:mPad(1)+n,mPad(2)+1:mPad(2)+m) = img;
    else
        % Do not pad
        img_pad = img;
    end
else
    n = numel(img);
    if(max(mPad > 0))
        
        % Pad image
        img_pad = zeros(2*mPad+n,1);
        img_pad(mPad+1:mPad+n) = img;
        img_pad((n+mPad+1):end) = fliplr(img(mPad:end));
        img_pad(1:mPad) = fliplr(img(1:mPad));
        
    else
        % Do not pad
        img_pad = img;
    end  
 
end

end
