function y = multiscale_problem
%% Construct 1D test problem Box and Gaussian
T = 24;
N = 64;

gausPos = 32;

y = zeros(1,N,1,T);
y(1,:,1,1) = gaussian_basis_1D(N,gausPos,20);
yy = decimate(y(1,:,1,1),4);
y(1,:,1,1:6) = repmat(yy(1,:,1)',[1,6]);
y(1,:,1,7:12) = repmat(yy(1,:,2)',[1,6]);
y(1,:,1,13:18) = repmat(yy(1,:,3)',[1,6]);
y(1,:,1,19:24) = repmat(yy(1,:,4)',[1,6]);

end
