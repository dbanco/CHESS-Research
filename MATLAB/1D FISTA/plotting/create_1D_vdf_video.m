%% Create VDF video
vdf_filename = 'overlap_vdf.gif';
x = sqrt(P.var_theta);
h = figure(6);
all_peaks1 = [];
all_peaks2 = [];
all_y1 = [];
all_y2 = [];
for i = 1:91
    im_vdf = squeeze(ring_data{1}.var_signal(i,:));
    
    % Find peaks in VDF
    peaks = find2peaks(im_vdf);
    if(numel(peaks) > 1)
        all_peaks2 = [all_peaks2, peaks(1)];
        all_peaks1 = [all_peaks1, peaks(2)];
        all_y2 = [all_y2, im_vdf(peaks(1))];
        all_y1 = [all_y1, im_vdf(peaks(2))];
    elseif (numel(all_peaks2)> 0)
        all_peaks2 = [all_peaks2, peaks(1)];
        all_y2 = [all_y2, im_vdf(peaks(1))];
    else
        all_peaks1 = [all_peaks1, peaks(1)];
        all_y1 = [all_y1, im_vdf(peaks(1))];
    end
    plot(x,im_vdf);
    hold on
    if( (numel(peaks) == 1) & (numel(all_peaks2)==0) )
        plot(x(peaks(1)),im_vdf(peaks(1)),'og')
        plot(x(all_peaks1),all_y1,'-g')
    elseif(numel(peaks)>1)
        plot(x(peaks(2)),im_vdf(peaks(2)),'og')
        plot(x(peaks(1)),im_vdf(peaks(1)),'or')
        plot(x(all_peaks1),all_y1,'-g')
        plot(x(all_peaks2),all_y2,'-r')
    else
        plot(x(peaks(1)),im_vdf(peaks(1)),'or')
        plot(x(all_peaks1),all_y1,'-g')
        plot(x(all_peaks2),all_y2,'-r')
    end
 
%     plot(x(peaks(2)),im_vdf(peaks(2)),'or')
    ylim([0,0.2])
    legend(num2str(ring_data{1}.az_spread(i)))
    drawnow
    pause(0.1)
    frameg = getframe(h);
    im = frame2im(frameg);

    [imind,cm] = rgb2ind(im,256,'nodither'); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,vdf_filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,vdf_filename,'gif','WriteMode','append','DelayTime',0.1); 
    end 
    cla
end