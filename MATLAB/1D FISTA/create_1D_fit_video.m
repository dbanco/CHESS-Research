%% Create VDF video
vdf_filename = 'overlap_fit.gif';
x = sqrt(P.var_theta);
h = figure(6);

for i = 1:91
    j = 1;
    fprintf('Image %i\n',i)
    fName = sprintf('fista_fit_%i_%i.mat',j,i);
    load(fullfile(fDir,fName))
    A0ft_stack = unshifted_basis_vector_ft_stack_norm2(P);
    img_fit = Ax_ft_1D(A0ft_stack,x_hat);
    plot(polar_vector)
    hold on
    plot(img_fit*norm(polar_vector))
    legend(num2str(ring_data{1}.az_spread(i)))
    ylim([0,0.5])
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