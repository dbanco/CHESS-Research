function y = loadMMPAD2D(ring,interp,topDir)
% Load MMPAD subset
if nargin < 3
    topDir = ['E:\Data\mmpad\ring',num2str(ring),'_zero'];
end
gif_filename = ['mmpad_ring',num2str(ring),'.gif'];
fName = 'mmpad_img_%i.mat';
% Load first image
f=load(fullfile(topDir,sprintf(fName,1)));
[theta,eta] = size(f.polar_image);

T = 546;
sub_theta = 1:theta;
sub_eta = 1:eta;
sub_T = 120;
t_step = 1;

h = figure(6);
h.Position = [-2.3486e+03 -318.6000 930 100.4000];

% Load subset of data
y = zeros(numel(sub_theta),numel(sub_eta),sub_T);
j = 1;
for i = t_step:t_step:sub_T
    f=load(fullfile(topDir,sprintf(fName,i)));
    
    x = f.polar_image(sub_theta,sub_eta);

    % Interpolate missing region
    if interp
        x(:,129:134) = (1:6)./7.*(x(:,135)-x(:,128)) + x(:,128);
    end
    y(:,:,j) = x;%/norm(x(:));

    j = j+1;

    h.Color = 'w';
    h.MenuBar = 'none';
    h.ToolBar = 'none';

    imagesc(x)
    axis off;

    drawnow
    pause(0.05)
    frame = getframe(gca);
    im = frame2im(frame);

    [imind,cm] = rgb2ind(im,256,'nodither'); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,gif_filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,gif_filename,'gif','WriteMode','append','DelayTime',0.1); 
    end 
    cla

end
end

