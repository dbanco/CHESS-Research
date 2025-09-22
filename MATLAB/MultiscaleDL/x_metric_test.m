% X-metric test script

K = 2;
J = 8;
N = 10;
T = 3;
X = zeros(1,N,K*J,T);
Xtrue = zeros(1,N,K*J,T);

% K,
Xtrue(1,5,5,:)=1;
X(1,6,5,:)=0.8;

Xtrue(1,7,13,:)=1;
X(1,5,14,:)=1;


dist = compute_x_metric(Xtrue,X,K,J);