function [I,stats] = cov_to_mvmi(V,x,y,m,N)

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

rstats = nargin > 3;

if rstats;
	if nargin < 5 || isempty(N), N = 1; end
	nx = length(x);
	ny = length(y);
	M = N*m;
	d  = nx*ny;
    Xm = d/M;           % chi^2 mean
	stats = struct('LR',struct('tstat',nan,'pval',nan));
	stats.LR.db = Xm;
end

z = 1:n; z([x y]) = [];
xz = [x,z];
yz = [y,z];

I = logdet(V(xz,xz)) + logdet(V(yz,yz)) - logdet(V(z,z)) - logdet(V);

if rstats
	stats.LR.tstat = I; % likelihood-ratio test statistic
    stats.LR.pval = 1-chi2cdf(M*stats.LR.tstat,d);
end
