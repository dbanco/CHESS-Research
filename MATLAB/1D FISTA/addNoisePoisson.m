function yn = addNoisePoisson(y,c)
%addNoisePoisson Summary of this function goes here
%   Detailed explanation goes here
yn = y + c*randn(size(y)).*sqrt(y);

end

