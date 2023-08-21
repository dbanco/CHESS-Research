function [u,v] = roundVelocity(u,v)

ratio = abs(u)./abs(v);

v(ratio>2) = 0;
u(ratio<0.5) = 0;

diagonals = (ratio < 2) & (ratio > 0.5);

magnitudes = max(u(diagonals),v(diagonals));

u(diagonals) = magnitudes.*sign(u(diagonals));
v(diagonals) = magnitudes.*sign(v(diagonals));



end