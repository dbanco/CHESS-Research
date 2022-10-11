function [c1min,c2min] = fareyApprox(input,lim)
% Approximate nearest rational number

a1 = 0;
a2 = 1;
b1 = 1;
b2 = 0;
c1 = a1;
c2 = a2;
err = abs(input-c1/c2);
c1min = c1;
c2min = c2;

while 1
    c1 = a1+b1;
    c2 = a2+b2;
    newErr = abs(input-c1/c2);
    if input*c2 > c1
        if lim < c2
            c1 = b1;
            c2 = b2;
            if newErr < err
                c1min = c1;
                c2min = c2;
            end
            return
        end
        a1 = c1;
        a2 = c2;
    elseif input*c2 == c1
        if lim >= c2
            if newErr < err
                c1min = c1;
                c2min = c2;
            end
            return
        end
        if b1 < b2
            c1 = a1;
            c2 = a2;
            if newErr < err
                c1min = c1;
                c2min = c2;
            end
            return
        end
        c1 = b1;
        c2 = b2;   
        if newErr < err
            c1min = c1;
            c2min = c2;
        end
        return
    else
        if lim < c2
            c1 = a1;
            c2 = a2;
            if newErr < err
                c1min = c1;
                c2min = c2;
            end
            return
        end
        b1 = c1;
        b2 = c2;
    end
    if newErr < err
        c1min = c1;
        c2min = c2;
        err = newErr;
    end
end

        
            
            
end
        
