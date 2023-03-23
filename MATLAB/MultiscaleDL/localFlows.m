function V = localFlows(Fy,Fx,Ft,win)
[Ny,Nx,Nt] = size(Fy);

V = zeros(floor([2,([Ny,Nx,Nt])/win]));
for t1 = 1:win:(Nt-win)
    for x1 = 1:win:(Nx-win)
        for y1 = 1:win:(Ny-win)
            % Index/pick out individual windows
            A = [vec(Fy(y1:(y1+win-1),x1:(x1+win-1),t1:(t1+win-1))),...
                 vec(Fx(y1:(y1+win-1),x1:(x1+win-1),t1:(t1+win-1)))];
            b = -vec(Ft(y1:(y1+win-1),x1:(x1+win-1),t1:(t1+win-1)) );
            % Solve linear system
            V( :, (y1+win-1)/win, (x1+win-1)/win, (t1+win-1)/win ) = A\b;
        end
    end
end

end