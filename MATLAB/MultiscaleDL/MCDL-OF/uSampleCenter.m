function du = uSampleCenter(d,m,center)
%decimatedDictionary Summary of this function goes here
%   Detailed explanation goes here

if m == 1
    du = d;
else
    [N1,N2] = size(d);
    if N2 > N1
        du = zeros(N1,N2*m);
        inds = [fliplr(center-m:-m:1),center:m:m*N2];
        du(:,inds) = d;
    else
        du = zeros(N1*m,N2);
        inds = [fliplr(center-m:-m:1),center:m:m*N1];
        du(inds,:) = d;
    end
%     norm(du(:))
%     figure(12)
%     hold on
%     plot(du(:,:,1))
    du = lowpassM(du,m,0);
    
end


end

