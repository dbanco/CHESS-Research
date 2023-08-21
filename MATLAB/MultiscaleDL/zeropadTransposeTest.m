
N = 4;
M = 3;
Total = N + M - 1;
a = [3;2;1];
b = padarray(a,N-1,0,'post');
x = [1 0 2 2 0 1]';

A = [1 2 3 0 0 0;
     0 1 2 3 0 0;
     0 0 1 2 3 0;
     0 0 0 1 2 3];

B = [3 0 0 0 1 2;
     2 3 0 0 0 1;
     1 2 3 0 0 0;
     0 1 2 3 0 0;
     0 0 1 2 3 0;
     0 0 0 1 2 3];

y1 = A*x;
y2 = B*x;
y2unpad = unpad(y2,2,'pre');
y3 = unpad(ifft(fft(b).*fft(x)),2,'pre');

y2unpadpad = padarray(y2unpad,2,0,'pre');
y3pad = padarray(y3,2,0,'pre');

Aty1 = A'*y1;
Bty2 = B'*y2;
Bty2unpad = B'*y2unpadpad;
Bty3unpad = ifft( conj(fft(b)).*fft(y3pad) );
% Bty3unpad = unpad( ifft( conj(fft(b)).*fft(y3pad) ),2,'post');

figure
subplot(2,1,1)
hold on
plot(y1,'o-')
plot(y2,'x-')
plot(y2unpad,'s-')
plot(y3,'v-')

subplot(2,1,2)
hold on 
plot(Aty1,'o-')
% plot(Bty2,'x-')
plot(Bty2unpad,'s-')
plot(Bty3unpad,'v-')

function x = unpad(xpad,M,shape)
switch shape
    case 'both'
        x = xpad(1+M:end-M,:);
    case 'post'
        x = xpad(1:end-M,:);
    case 'pre'
        x = xpad(1+M:end,:);
end
end