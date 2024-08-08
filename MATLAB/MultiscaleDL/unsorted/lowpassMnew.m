function y = lowpassM(x,m,trans)
N = numel(x);
M = 2*m-1;
fout = zeros(1,N);
% Convolution kernel (linear interpolation)
if m == 1
    y = x;
else
    f = zeros(1,M);
    f(1:m) = linspace(1/m,1,m);
    f(m+1:end) = f(m-1:-1:1);
    f = f/norm(f);
    fout(1) = f(m);
    fout(2:m) = f(m-1:-1:1);
    fout(N:-1:(N-m+2)) = f(1:m-1);
    y = ifft( fft(x).*fft(fout));
end




end