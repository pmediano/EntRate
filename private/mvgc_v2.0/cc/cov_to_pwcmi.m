function [I,stats] = cov_to_pwcmi(V,m,N)

% pairwise-conditional MIs (partial correlation!)

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

rstats = nargin > 1;

if rstats;
	if nargin < 3 || isempty(N), N = 1; end
	M = N*m;
    d  = 1;   % chi^2 df and F df1
    Xm = d/M; % chi^2 mean
	stats = struct('LR',struct('tstat',nan(n),'pval',nan(n)));
	stats.LR.db = Xm;
end

LDV = logdet(V);

LDVI = zeros(n,1);
for i = 1:n
    oi = 1:n; oi(i) = [];
    LDVI(i) = logdet(V(oi,oi));
end

I = nan(n);

for i = 1:n
	for j = i+1:n
        oij = 1:n; oij([i j]) = [];
        I(i,j) = LDVI(i) + LDVI(j) - logdet(V(oij,oij)) - LDV;
        I(j,i) = I(i,j);
	end
end

if rstats
	stats.LR.tstat = I; % likelihood-ratio test statistic
    stats.LR.pval  = 1-chi2cdf(M*I,d);
end
