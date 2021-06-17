function mre = maxrelerr(x,y)

assert(isequal(size(x),size(y)),'Arrays must have the same dimensions');

d = abs(y-x);
x = abs(x);
x(x <= 2*eps) = 1; % minimum detectable difference between x and a value close to x is O(x)*eps.
relerr = d./x;
mre = max(relerr(:));
