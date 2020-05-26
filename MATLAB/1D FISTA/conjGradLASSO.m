function xk = conjGradLASSO( A0ft_stack, xk, yk, vk, b, zMask, rho, maxGradIters, gradTol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

bnorm = norm(b(:));
r10 = AtR_ft_1D(A0ft_stack,b - forceMaskToZero(Ax_ft_1D(A0ft_stack,xk), zMask));
r20 = rho*((yk - vk) - xk);
r1k = r10;
r2k = r20;
p10 = r10;
p1k = p10;
p20 = r20;
p2k = p20;
for k = 1:maxGradIters
    alpha1k = sum(sum(r1k.*r1k))/sum(sum( p1k.*AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,p1k), zMask)) ));
    alpha2k = sum(sum(r2k.*r2k))/sum(sum(p2k.*p2k));
    xk = xk + alpha1k*p1k + alpha2k*p2k;
    r1kp1 = r1k - alpha1k*AtR_ft_1D(A0ft_stack,forceMaskToZero(Ax_ft_1D(A0ft_stack,p1k), zMask));
    r2kp1 = r2k - alpha2k*p2k;
    
    if (norm(r1kp1) + norm(r2kp1)) < gradTol 
       break; 
    end
    
    beta1k = sum(sum(r1kp1.*r1kp1))/sum(sum(r1k.*r1k));
    beta2k = sum(sum(r2kp1.*r2kp1))/sum(sum(r2k.*r2k));
    p1k = r1kp1 + beta1k*p1k;
    p2k = r2kp1 + beta2k*p2k;
    r1k = r1kp1;
    r2k = r2kp1;
end

end

