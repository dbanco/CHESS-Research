function out = centralDiff(in,dim)
%centralDiff Computes central difference along dim
N = size(in,dim);

in_n = unfold(in,dim);

out = zeros(size(in_n));
out(1,:) = in_n(2,:) - in_n(1,:);
out(N,:) = in_n(N,:) - in_n(N-1,:);
out(2:(N-1),:) = (in_n(3:N,:) - in_n(1:N-2,:))/2;

out = fold(out,size(in),dim);
end