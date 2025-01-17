function y = tooth(N,width)

y = zeros(N,1);
center = (N+1)/2;

if mod(width,2)
    y(center-floor(width/2): center+floor(width/2)) = linspace(0,1,width)/norm(linspace(0,1,width));
else
    y(center-floor(width/2): center+floor(width/2)-1) = linspace(0,1,width)/norm(linspace(0,1,width));
end

end