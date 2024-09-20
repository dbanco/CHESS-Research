% Detector setup

center = [1000,1000];

twotheta = 7.61;
radius = 610*tan(pi*7.61/180)/0.0748 %mm 

% Ti-7
a = 2.92255e-8;
b = 2.932e-8;
c = 4.67133e-8;
detectDist = 610;
indices = [3,1,1; 2,2,2; 4,0,0];

% c103
a = 0.33139e-10;
b = 0.33139e-10;
c = 0.33139e-10;
detectDist = 884.3593298614039;
indices = [1,1,0; 1,1,2; 1,2,3];

%  constants
eV = 1.60218e-19;
energy = 41.1*eV;
sol = physconst('LightSpeed');
h = 6.62607015e-34;

lam = h*sol/energy;

N = size(indices,1);
thetas = zeros(N,1);

for i = 1:size(indices,1)
    d = sqrt(indices(i,1)^2 + indices(i,2)^2 + indices(i,3)^2);
    thetas(i) = asin(lam/a*d)/2;
    rad(i) = detectDist*tan(thetas(i));
end



