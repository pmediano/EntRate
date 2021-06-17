function [I,stats] = cov_to_gwcmi(V,group,m,N)

% For each pair of groups, returns mutual information
% between them, conditional on rest of system.

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

g = check_group(group,n);

rstats = nargin > 2;

if rstats;
	if nargin < 4 || isempty(N), N = 1; end
	M = N*m;
    d  = nan(g);  % chi^2 df and F df1
    for a = 1:g
		na = length(group{a});
		for b = 1:g
			if b == a, continue; end
			nb = length(group{b});
			d(a,b)  = na*nb;
		end
	end
    Xm = d/M;           % chi^2 mean
	stats = struct('LR',struct('tstat',nan(g),'pval',nan(g)));
	stats.LR.db = Xm;
end

LDV = logdet(V);

LDVG = zeros(g,1);
for a = 1:g
    goa = 1:n; goa(group{a}) = [];
    LDVG(a) = logdet(V(goa,goa));
end

I = nan(g);

for a = 1:g
	for b = a+1:g
        goab = 1:n; goab([group{a} group{b}]) = [];
        I(a,b) = LDVG(a) + LDVG(b) - logdet(V(goab,goab)) - LDV;
        I(b,a) = I(a,b);
	end
end

if rstats
	stats.LR.tstat = I; % likelihood-ratio test statistic
    stats.LR.pval  = 1-chi2cdf(M*I,d);
end
