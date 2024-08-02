[yn,y,K,J,N,M,T,Xtrue,Dtrue,scales] = gaus_linear_osc_signal_matched_small_zpad2_center(0);

close all

% figure
% hold on
% plot(Dtrue(:,:,1))
% plot(Dtrue(:,:,2))
 
AD = reSampleCustomArray(N,Dtrue,scales);
AtAD = reSampleTransCustomArray(M,AD,scales);

% Regular dictionary
figure(1)
hold on
plot(AtAD(:,:,1))
plot(AtAD(:,:,2))

center = (M+1)/2;
ADcenter = reSampleCustomArrayCenter(N,Dtrue,scales,center);
AtADcenter = reSampleTransCustomArrayCenter(M,ADcenter,scales,center);

figure(2)
hold on
title('Centered Dict Atom 1')
for i = 1:J
    plot(ADcenter(:,:,i))
end
plot([center,center],[0,0.6],'o-')

figure(3)
hold on
title('Centered Dict Atom 2')
for i = (J+(1:J))
    plot(ADcenter(:,:,i)/max(ADcenter(:,:,i)))
end
plot([center,center],[0,0.6],'o-')

figure(4)
title('AtAD 1 and 2')
hold on
plot(AtADcenter(:,:,1))
plot(AtADcenter(:,:,2))

D0 = zeros(size(D));
D0(:,23,:) = 1;
AD = reSampleCustomArray(N,D0,scales);
AtAD = reSampleTransCustomArray(M,AD,scales);
ADcenter = padarray(reSampleCustomArrayCenter(N,D0,scales,center),[0,M-1,0,0],0,'post');
AtADcenter = reSampleTransCustomArrayCenter(M,ADcenter,scales,center);

figure(5)
title('AtAD Regular and Centered (Atom 1)')
hold on
plot(AtAD(:,:,1),'o-')
plot(AtADcenter(:,:,1),'-x')
