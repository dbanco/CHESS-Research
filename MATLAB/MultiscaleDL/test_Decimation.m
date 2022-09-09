clear all
close all
% Construct circulant X,x
x1 = zeros(64,1);
x1(10,1) = 1;
X1 = zeros(64);
for i = 1:64
    X1(:,i) = circshift(x1,i);
end

x2 = zeros(64,1);
x2(40,1) = 1;
X2 = zeros(64);
for i = 1:64
    X2(:,i) = circshift(x2,i);
end

% Construct circulant A,a
a = zeros(64,1);
a(1:10) = 1;
a(11:20) = 1.5;
A = zeros(64);
for i = 0:63
    A(:,i+1) = circshift(a,i);
end

[Da,D] = decimatedDictionary(a,3);

% plotDictionary([X1 X2]*Da)
% title('decimate a')
% plotDictionary([X*D*a);
% title('decimate x')

%% Arbitrary circulant matrix
D1 = D(1:64,:);
% D1(1:2:63,:) = 0;
D2 = D(65:128,:);
D3 = D(129:192,:);
c = rand(64,1);
C = zeros(64);
for i = 1:64
    C(:,i) = circshift(c,i);
end
DFT = dftmtx(64);
figure(1)
imagesc(D1*C)
figure(11)
imagesc(D2*C)
figure(111)
imagesc(D3*C)
figure(2)
imagesc(abs(DFT*C*D1*DFT'))
figure(3)
imagesc(abs(DFT*C*D2*DFT'))
figure(4)
imagesc(abs(DFT*C*D3*DFT'))