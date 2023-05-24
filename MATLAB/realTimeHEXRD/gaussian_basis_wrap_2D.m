function b = gaussian_basis_wrap_2D( N_y,mean_y,std_y,...
                                     N_x,mean_x,std_x,scaling)
%gaussian_basis_wrap_2D Generates 2d gaussian basis function
%
% Inputs:
% N_y- image size in radial direction
% mean_y - radial mean of gaussian basis function
% std_y - radial variance of gaussian basis function
% N_x - image size in theta direction
% mean_x - mean of gaussian basis function in theta
% std_x - standard deviation of gaussian basis function in theta
%
% Outputs:
% b - (N_y x N_x) array

wrapN = @(x, N) (1 + mod(x-1, N));

% Compute y distances with wrapping at boundaries
idx = 1:N_y;
mean_y = wrapN(mean_y,N_y);
upDiff = ceil(mean_y) - mean_y;
downDiff = mean_y - floor(mean_y);
upDist = idx - 1 + upDiff;
downDist = flip( idx - 1 + downDiff);
disty1 = circshift(upDist,floor(mean_y));
disty2 = circshift(downDist,floor(mean_y));
if mean_y == floor(mean_y); disty1 = circshift(disty1,-1); end

% Compute x distances with wrapping at boundaries
idx = 1:N_x;
mean_x = wrapN(mean_x,N_x);
upDiff = ceil(mean_x) - mean_x;
downDiff = mean_x - floor(mean_x);
upDist = idx - 1 + upDiff;
downDist = flip( idx - 1 + downDiff);
distx1 = circshift(upDist,floor(mean_x));
distx2 = circshift(downDist,floor(mean_x));
if mean_x == floor(mean_x); distx1 = circshift(distx1,-1); end

b1 = zeros(N_y,N_x);
b2 = zeros(N_y,N_x);
for i = 1:N_y
    for j = 1:N_x
        b1(i,j) = exp(-disty1(i).^2/(2*std_y^2) -...
                       distx1(j).^2/(2*std_x^2));
        b2(i,j) = exp(-disty2(i).^2/(2*std_y^2) -...
                       distx2(j).^2/(2*std_x^2));
        b3(i,j) = exp(-disty1(i).^2/(2*std_y^2) -...
                       distx2(j).^2/(2*std_x^2))';
        b4(i,j) = exp(-disty2(i).^2/(2*std_y^2) -...
                       distx1(j).^2/(2*std_x^2))';
                   
    end
end
b = b1+b2+b3+b4;

if floor(mean_x) == mean_x
    b(:,mean_x) = b3(:,mean_x) + b4(:,mean_x);
end
if floor(mean_y) == mean_y
    b(mean_y,:) = b3(mean_y,:) + b4(mean_y,:);
end
if floor(mean_x) == mean_x && floor(mean_y) == mean_y
    b(mean_y,mean_x) = b1(mean_y,mean_x);
end

% Not going to account for signal that wraps more than once
% for k = 1:floor(6*std_x/N_x)
%     for i = 1:N_y
%         for j = 1:N_x
%             b1(i,j) = exp(-(disty1(i)+k*N_y).^2/(2*std_y^2) -...
%                            (distx1(j).^2/(2*std_x^2))';
%             b2(i,j) = exp(-disty2(i).^2/(2*std_y^2) -...
%                            distx2(j).^2/(2*std_x^2))';
%         end
%     end
%     b = b + b1 + b2;
% end

if nargin > 6
    switch scaling
        case '2-norm'
            b = b/norm(b(:));
        case '1-norm'
            b = b/sum(abs(b(:)));
        case 'max'
            b = b/max(b(:));
        case 'rms'
            b = b/sqrt( sum(b(:).^2)/N_x );
        otherwise
            b = b/(sqrt(std_x*std_y)*2*pi);
    end
else
    b = b/(sqrt(std_x*std_y)*2*pi);
end

end

