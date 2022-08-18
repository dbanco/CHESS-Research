function [ peaks ] = find2peaks( y )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
center = y(4:end-3);
diff1 = center - y(1:end-6);
diff2 = center - y(2:end-5);
diff3 = center - y(3:end-4);

diff4 = center - y(5:end-2);
diff5 = center - y(6:end-1);
diff6 = center - y(7:end);

peaks = find(((diff1 > 0) &...
              (diff2 > 0) &...
              (diff3 > 0) &...
              (diff4 > 0) &...
              (diff5 > 0) &...
              (diff6 > 0) ), 2) + 3;
end

