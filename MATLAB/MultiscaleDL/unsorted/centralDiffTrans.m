function out = centralDiffTrans(in,dim)
%centralDiffTrans Computes transposed central difference along dim
N = size(in,dim);
out = zeros(N);
in_n = unfold(in,dim);

out(1,:) = -in_n(2,:) - in_n(1,:);
out(N,:) = in_n(N,:) + in_n(N-1,:);
out(2:(N-1)) = (-in_n(3:N,:) + in_n(1:N-2,:))/2;
end