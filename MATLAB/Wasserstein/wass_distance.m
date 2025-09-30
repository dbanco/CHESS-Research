function wass_dist = wass_distance(Xtrue,X,K,J)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Drop singleton first dim
Xtrue = squeeze(Xtrue);  % [N x (K*J) x T]
X  = squeeze(X);   % [N x (K*J) x T]

N = size(Xtrue,1);       % number of "shifts" (M should equal N?)
T = size(Xtrue,3);

% Reshape into [K x J x N x T]
tmp1 = zeros(K,J,N,T);
tmp2 = zeros(K,J,N,T);
for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    tmp1(k,:,:,:) = permute(Xtrue(:,i1:i2,:),[2,1,3]);
    tmp2(k,:,:,:) = permute(X(:,i1:i2,:),[2,1,3]);
end
Xtrue = tmp1;
X  = tmp2;

Xtrue = squeeze(sum(Xtrue,1));
X = squeeze(sum(X,1));

D = grid_distance_matrix(J,N);

wass_dist = 0;
for t = 1:T
    % (num2str(t))
    vec1 = Xtrue(:,:,t);
    vec2 = X(:,:,t);

    % figure
    % subplot(2,1,1)
    % imagesc(Xtrue(:,:,t))
    % subplot(2,1,2)
    % imagesc(X(:,:,t))

    wass_dist = wass_dist + sinkhornKnoppTransport(vec1(:),vec2(:), 10, D);
end