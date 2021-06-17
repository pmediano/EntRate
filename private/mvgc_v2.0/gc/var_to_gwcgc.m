function [F,stats] = var_to_gwcgc(A,V,group,whichstats,X,regmode)

% inter-group (conditional) GCs

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(group,n);

stats1 = false;
stats2 = false;
if nargin > 3 && ~isempty(whichstats)
	switch lower(whichstats)
		case 'single', stats1 = true;
		case 'dual',   stats2 = true;
		case 'both',   stats1 = true; stats2 = true;
		case 'none'    % do nothing
		otherwise error('bad stats specification');
	end
end

stats = [];
if stats1 | stats2;
	if     stats2
		assert(nargin > 4 && ~isempty(X) && ~isscalar(X),'must supply time series data for dual-regression stats');
		assert(nargin > 5,'must supply regression mode (same as for parameter estimates) for dual-regression stats');
        [~,m,N] = size(X);
	elseif stats1
		assert(nargin > 4 && ~isempty(X),'must supply number of observations and trials for single-regression stats');
		if isscalar(X)
			assert(nargin > 3 && ~isempty(regmode),'must supply number of observations and trials for single-regression stats');
			m = X;
			N = regmode;
			clear X regmode;
		else
			[~,m,N] = size(X);
		end
	end
    M  = N*(m-p); % chi^2 scaling factor = effective number of observations
    d  = nan(g);  % chi^2 df and F df1
    d2 = nan(g);  % F df2
    for b = 1:g
		ny = length(group{b});
		for a = 1:g
			if a == b, continue; end
			nx = length(group{a});
			d(a,b)  = p*nx*ny;
			d2(a,b) = nx*(M-p*n)-1;
		end
		K  = d2./d;
		Xm = d/M;
		Fm = d./(d2-2);
	end
    statsf = struct('tstat',nan(g),'pval',nan(g));
    if stats1
		stats.single = struct('F',statsf,'LR',statsf);
		stats.single.F.db  = Fm;
		stats.single.LR.db = Xm;
	end
    if stats2
		stats.dual = struct('F',statsf,'LR',statsf);
		stats.dual.F.db  = Fm;
		stats.dual.LR.db = Xm;
	end
    clear statsf;
end

F = NaN(g);

for b = 1:g
	y = group{b};
    r = 1:n; r(y) = []; % omit group b

	[~,VR,rep] = var2riss(A,V,y,r);
    if sserror(rep,b), continue; end % check DARE report, bail out on error

	if stats2
		[~,VR2] = tsdata_to_var(X(r,:,:),p,regmode);
	end

    for a = 1:g
        if a == b, continue; end
        x = group{a};
        xr = findin(x,r); % indices of group{a} in r

        F(a,b) = logdet(VR(xr,xr))-logdet(V(x,x));

		if stats1
			stats.single.F.tstat(a,b)  = trace(VR(xr,xr))/trace(V(x,x)) - 1;  % F-test statistic
			stats.single.LR.tstat(a,b) = F(a,b);                              % likelihood-ratio test statistic
		end

		if stats2
			stats.dual.F.tstat(a,b)  = trace(VR2(xr,xr))/trace(V(x,x)) - 1;   % F-test statistic
			stats.dual.LR.tstat(a,b) = logdet(VR2(xr,xr)) - logdet(V(x,x));   % likelihood-ratio test statistic
		end
    end
end

if stats1
    stats.single.F.pval  = 1-fcdf(K.*stats.single.F.tstat,d,d2);
    stats.single.LR.pval = 1-chi2cdf(M*stats.single.LR.tstat,d);
end

if stats2
    stats.dual.F.pval    = 1-fcdf(K.*stats.dual.F.tstat,d,d2);
    stats.dual.LR.pval   = 1-chi2cdf(M*stats.dual.LR.tstat,d);
end
