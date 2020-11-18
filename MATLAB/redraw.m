function redraw(frame)
load(['E:\PureTiRD_full\ring3_zero\mmpad_img_',num2str(frame),'.mat'])
imagesc(polar_image)
title(sprintf('%i',frame))
end