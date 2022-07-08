function f = plotMultiscaleDictionary(D,Phi)
Df = fft(D);
Phif = fft2(Phi,size(Phi,1),size(Phi,2));
PhifGf = scaleDict(Phif,Df);
PhiG = real(ifft2(PhifGf,'symmetric'));

K1 = size(D,3);
K2 = size(Phi,3);

f = figure;
rows = K1;
cols = K2;
k = 1;
for i = 1:K1
    for j = 1:K2
        subplot(rows,cols,j+K2*(i-1))
        plot(PhiG(:,:,k)) 
        set(gca,'YTickLabel',[]);
    %     set(gca,'XTickLabel',[]);
    k = k + 1;
    end
end
end

