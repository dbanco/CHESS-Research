function plotVDF(X)

figure;
vdf = squeeze(sum(sum(X,1),2));
imagesc(vdf)

end

