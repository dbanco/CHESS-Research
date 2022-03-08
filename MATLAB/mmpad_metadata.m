% MMPAD metadata
mmpad_dims = [396,265];
rings = {'020','112','021','004'};
d_space = [1.27773,1.24845,1.23272,1.17138];
two_theta = [6.89132,7.05316,7.14328,7.5179];

pixel_side = 150e-6;
gap_width = 0.75e-3;

mmpad_dims*pixel_side

detec_dist = 4.6;
detect_angle = 14.6;

% Compute circumeference of rings
radius = 4.6*tan(pi/90*two_theta/2);
circum = 2*pi*radius;

% Assume pixels approximately lie on circumeference
pixel_angle = pixel_side./circum*360;

% Compute 2theta angle
pixel_delta = 1:30;


function rad_angle = pixel2angle_radius(pixel_delta,init_radius,detec_dist,init_angle)
    rad_angle = arctan((init_radius + pixel_delta)/detec_dist) - init_angle;
end