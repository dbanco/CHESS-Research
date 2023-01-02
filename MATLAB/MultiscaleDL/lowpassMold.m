function y = lowpassMold(x,m,trans)
N = numel(x);
M = 2*m-1;

% Convolution kernel (linear interpolation)
if m == 1
    f = 1;
else
    f = zeros(1,M);
    f(1:m) = linspace(1/m,1,m);
    f(m+1:end) = f(m-1:-1:1);
end
f = f/norm(f);

fout = zeros(1,N);
if m > 1
    if trans 
        fout(1) = f(1);
        fout(N-(M-2):N) = f(M:-1:2);
    else
        fout(1:M) = f;
    end
    y = ifft( fft(x).*fft(fout));
else
    y = x;
end

end