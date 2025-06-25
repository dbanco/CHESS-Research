function du = uSample(d,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

if m == 1
    du = d;
else
    [N1,N2] = size(d);
    if N2 > N1
        du = zeros(N1,N2*m);
        du(:,m*(1:N2)) = d;
    else
        du = zeros(N1*m,N2);
        du(m*(1:N1),:) = d;
    end

    du = lowpassM(du,m,0);
    
end


end

