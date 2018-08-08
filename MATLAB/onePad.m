function img_pad = onePad( img, oPad )
%zeroPad Zero pads image

[n,m] = size(img);

if(max(oPad > 0))
    % Pad image
    img_pad = ones(2*oPad(1)+n,2*oPad(2)+m);
    img_pad(oPad(1)+1:oPad(1)+n,oPad(2)+1:oPad(2)+m) = img;
else
    % Do not pad
    img_pad = img;
end

end

