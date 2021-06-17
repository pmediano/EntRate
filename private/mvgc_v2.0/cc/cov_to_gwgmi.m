function [I,stats] = cov_to_gwgmi(V,group,m,N)

% For each group, returns multi-information of group conditional
% on rest of system.

[n,n1] = size(V);
assert(n1 == n,'Covariance matrix must be square');

[g,gsiz] = check_group(group,n);

rstats = nargin > 2;

if rstats;
	if nargin < 4 || isempty(N), N = 1; end
	M = N*m;
    d  = nan(g,1);  % chi^2 df and F df1
    for a = 1:g
		na = length(group{a});
		d(a)  = (na*(na-1))/2;
	end
    Xm = d/M;           % chi^2 mean
	stats = struct('LR',struct('tstat',nan(g,1),'pval',nan(g,1)));
	stats.LR.db = Xm;
end

I = nan(g,1);

LDV  = logdet(V);
for a = 1:g
    goa = 1:n; goa(group{a}) = [];
    I(a) = -LDV - (gsiz(a)-1)*logdet(V(goa,goa));
    for i = group{a}
        igoa = [i goa];
        I(a) = I(a) + logdet(V(igoa,igoa));
    end
end

if rstats
	stats.LR.tstat = I; % likelihood-ratio test statistic
    stats.LR.pval  = 1-chi2cdf(M*I,d);
end
