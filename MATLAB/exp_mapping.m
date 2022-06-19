% exponential map
x = 1; y = 2; lambda = 3;
v = [x;y;lambda];

s = exp(lambda);
X = (1-1/s)/lambda; Y = (1/s-1+lambda)/lambda^2;
V = [X 0;
     0 X];
algv = [ 0 0 x;...
         0 0 y;...
         0 0 -lambda];
  
Ghat = [ 1 0 X*x;...
         0 1 X*y;...
         0 0 1/s]
expG = expm(algv)

% logGhat = log(Ghat)
% logexpG = log(expG)
% 
% algv2 = [ 0 -1 x;...
%          1 0  y;...
%          0 0 -lambda];
% exp(x)