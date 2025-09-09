function dist = compute_x_metric(Xtrue, X, K, J)
% compute_x_metric
%   Compute distances between active coefficients in Yhat and Xtrue
%   across time steps, partitioning scale/shift space.
%
% Inputs:
%   Xtrue : [1 x N x (K*J) x T] ground truth coefficients (binary)
%   Yhat  : [1 x N x (K*J) x T] recovered coefficients
%   K     : #atoms
%   J     : #scales
%   M     : #shifts
%   T     : #time steps
%
% Output:
%   dists : {T x K} cell array
%           Each cell contains vector of distances for that (time, atom).

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

    truthCoords = zeros(2,K,T); %scale, shift
    distMaps = zeros(K,J,N,T);
    maskMaps = zeros(K,J,N,T);
    % 1: Identify truth coordinates and compute distance maps
    for t = 1:T
        for k = 1:K
            Xslice = squeeze(Xtrue(k,:,:,t));
            [jTrue, sTrue] = find(Xslice > 0.5);
            truthCoords(:,k,t) = [jTrue, sTrue];
            distMaps(k,:,:,t) = distance_map(jTrue,sTrue,J,N);
        end
    end

    % Mask distance maps
    minDistMap = min(distMaps,[],1);

    % Compute metric for each k
    Xsum = sum(X,1);
    dist = sum(minDistMap.*Xsum,'all');
end
