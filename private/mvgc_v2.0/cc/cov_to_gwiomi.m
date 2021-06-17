function [I,stats] = cov_to_gwiomi(V,group,m,N)

% For each group, returns MI between group and rest of system

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

g = check_group(group,n);

rstats = nargin > 2;

if rstats;
	if nargin < 4 || isempty(N), N = 1; end
	M = N*m;
	d  = nan(g,1);  % chi^2 df and F df1
	for a = 1:g
		na = length(group{a});
		nb = n-na;
		d(a)  = na*nb;
	end
	Xm = d/M;           % chi^2 mean
	stats = struct('LR',struct('tstat',nan(g,1),'pval',nan(g,1)));
	stats.LR.db = Xm;
end

LDV = logdet(V);

I = nan(g,1);
for a = 1:g
	ga = group{a};
	gb = 1:n; gb(ga) = [];
	I(a) = logdet(V(ga,ga)) + logdet(V(gb,gb)) - LDV;
end

if rstats
	stats.LR.tstat = I; % likelihood-ratio test statistic
    stats.LR.pval  = 1-chi2cdf(M*I,d);
end
