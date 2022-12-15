function du = uSample(d,m)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

if m == 1
    du = d;
else
    [N1,N2] = size(d);
    du = zeros(N1,N2*m);
    du(:,m*(1:N2)) = d;
%     norm(du(:))
%     figure(12)
%     hold on
%     plot(du(:,:,1))
    du = lowpassM(du,m,0);
    
end


end

