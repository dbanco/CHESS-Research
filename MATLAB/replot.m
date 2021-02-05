function replot(frame)
load(['E:\PureTiRD_full\ring4_zero\mmpad_img_',num2str(frame),'.mat'])
polar_vector = squeeze(sum(polar_image,1));
plot(polar_vector)
ylim([0 1e6])
text(numel(polar_vector),5e5,sprintf('%i',frame))
end