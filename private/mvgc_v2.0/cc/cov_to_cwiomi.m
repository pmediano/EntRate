function [I,stats] = cov_to_cwiomi(V,m,N)

% Returns MIs between each source and rest of system

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

rstats = nargin > 1;

if rstats;
	if nargin < 3 || isempty(N), N = 1; end
	M = N*m;
	d = n*(n-1);
	Xm = d/M;           % chi^2 mean
	stats = struct('LR',struct('tstat',nan(n,1),'pval',nan(n,1)));
	stats.LR.db = Xm;
end

LDV = logdet(V);

I = nan(n,1);
for i = 1:n
	oi = 1:n; oi(i) = [];
	I(i) = log(V(i,i)) + logdet(V(oi,oi)) - LDV;
end

if rstats
	stats.LR.tstat = I; % likelihood-ratio test statistic
    stats.LR.pval  = 1-chi2cdf(M*I,d);
end
