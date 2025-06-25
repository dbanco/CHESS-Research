% Center Penalty test script
M = 9;
c = (M - 1) / 2;                % Center of support
w = ((0:M-1) - c);          % Distance weights
w = w(:);                      % Ensure column vector
wwT = w*w';