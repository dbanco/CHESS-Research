function [ theta_project, theta_domain, theta_ind ] = azimuthal_projection( img,center,r,theta_domain )
%AZIMUTHAL_PROJECTION Summary of this function goes here
%     Interpolates along a line in the radial direction
%     
%     inputs
%             image         image
%             center      center of the rings in the image
%             r           radius of ring of interest     
%             num_theta   number of samples along theta
%             
%     outputs
%             theta_project  image values at points along line
%             theta_domain   domain over which image values are defined
%
    [n,m] = size(img);
    if(center==0) 
        center = [round(n/2.0), round(m/2.0)];
    end

    
    theta_project = zeros(size(theta_domain));
    
    for tidx = 1:numel(theta_domain)
        % identify surrounding four points
        x = r*cos(theta_domain(tidx)) + center(1);
        y = r*sin(theta_domain(tidx)) + center(2);    
        x1 = floor( x );
        x2 = ceil(  x );
        y1 = floor( y );
        y2 = ceil(  y );
  
        % make sure we are in image
        if( (y1 <= n) && (y2 <= n) && (x1 > 0) && (x2 > 0) &&...
            (x1 <= m) && (x2 <= m) && (y1 > 0) && (y2 > 0) )
            if((x2-x1 == 0) && (y2-y1 == 0))
                theta_project(tidx) = img(y1,x1);
            elseif((x2-x1) == 0)
                theta_project(tidx) = img(y1,x1) +...
                                (img(y2,x2)-img(y1,x1))*(y-y1)/(y2-y1);
            elseif((y2-y1) == 0)
                theta_project(tidx) = img(y1,x1) +...
                                (img(y2,x2)-img(y1,x1))*(x-x1)/(x2-x1);
            else
                
                % interpolate
                a = [x2-x; x-x1];
                Q = [img(y1,x1), img(y1,x2);
                     img(y2,x1), img(y2,x2)];     
                b = [y2-y; y-y1];
                theta_project(tidx) = a'*Q*b/((x2-x1)*(y2-y1));
            end
        else
            theta_project(tidx) = NaN;
        end
    end
    removeVals = isnan(theta_project);
    theta_project(removeVals) = [];
    theta_ind = ~removeVals;
end

