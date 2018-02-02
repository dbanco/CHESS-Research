%% Identify XRD Rings
detect_dist = 3293.09;
radii =   [732,   766,   890]; 
indices = [3,1,1; 2,2,2; 4,0,0];
lambda = 2.6379627e-10;
theta1 = asin(radii(1)/detect_dist)/2;
% Use 311 ring to compute lattice parameter
a = sqrt(11/4*lambda^2/sin(theta1/2)^2);

%% Compare measured angle of ring to theoretical angle given by miller indices
for i = 1:numel(radii)
    theta1 = asin(radii(i)/detect_dist)/...
             (2*pi)*360/2;
    theta2 = asin(sqrt(lambda^2/a^2*(indices(i,1)^2 + indices(i,2)^2 + indices(i,3)^2)))/...
             (2*pi)*360/2;
    fprintf('Image: %2.4f, %i%i%i: %2.4f \n',theta1,...
            indices(i,1),indices(i,2),indices(i,3),theta2) 
end

