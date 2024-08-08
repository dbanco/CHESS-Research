function out = velCon(x,K,u,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[N, KJ, T] = size(x);

epsilon = 1e-16;
Vtm1 = zeros(N,KJ,T-1);
Vt1 = zeros(N,KJ,T-1);
noV = zeros(N,KJ,T);

smoothness = 1;
maxIt = 1;

if nargin < 3
    [u,v,~,~,~]  = computeHornSchunkDict(x,K,smoothness,maxIt);
end

Vtm1(abs(u(:,:,2:T))>epsilon) = -1;
Vtm1(abs(v(:,:,2:T))>epsilon) = -1;

% Place 1 where velocities point towards
for i = 1:N
    for j = 1:KJ
        for t = 2:T
            % Check x-direction
            if u(i,j,t) > epsilon
                jj = j+1;
                if jj > T-1, jj = T-1; end
            elseif u(i,j,t) < -epsilon
                jj = j-1;
                if jj <1, jj = 1; end
            else
                jj = j;
            end
            % Check y-direction
            if v(i,j,t) > epsilon
                ii = i+1;
                if ii > T-1, ii = T-1; end
            elseif v(i,j,t) < -epsilon
                ii = i-1;
                if ii < 1, ii = 1; end
            else
                ii = i;
            end
            % Set value
            if ~((ii == i) && (jj == j))
                Vt1(ii,jj,t-1) = Vt1(ii,jj,t-1) + 1;
            end
        end
    end
end

% Define where velocitys are not to penalize more heavily
for t = 2:T
    [is,js] = find((Vtm1(:,:,t-1)>0));
    noV(is,js,t-1) = 1;
    [is,js] = find((Vt1(:,:,t-1)>0));
    noV(is,js,t) = 1;
end


% Compute velocity constraint
tsum = sum(sum( (Vt1.*x(:,:,2:T) + Vtm1.*x(:,:,1:T-1) ),1),2);
plot(squeeze(tsum))

out = norm(tsum(:)).^2;

end