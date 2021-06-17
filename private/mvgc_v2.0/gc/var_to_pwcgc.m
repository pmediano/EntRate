function [F,stats] = var_to_pwcgc(A,V,whichstats,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

stats1 = false;
stats2 = false;
if nargin > 2 && ~isempty(whichstats)
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
		assert(nargin > 3 && ~isempty(X) && ~isscalar(X),'must supply time series data for dual-regression stats');
		assert(nargin > 4,'must supply regression mode (same as for parameter estimates) for dual-regression stats');
        [~,m,N] = size(X);
	elseif stats1
		assert(nargin > 3 && ~isempty(X),'must supply number of observations and trials for single-regression stats');
		if isscalar(X)
			assert(nargin > 4 && ~isempty(regmode),'must supply number of observations and trials for single-regression stats');
			m = X;
			N = regmode;
			clear X regmode;
		else
			[~,m,N] = size(X);
		end
	end
    M  = N*(m-p);  % chi^2 scaling factor = effective number of observations
    d  = p;        % chi^2 df and F df1
    Xm = d/M;      % chi^2 mean
    d2 = M-p*n-1;  % F df2
    K  = d2/d;     % F scaling factor
    Fm = d/(d2-2); % F mean
    statsf = struct('tstat',nan(n),'pval',nan(n));
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

DV = diag(V);
LDV = log(DV);

F = nan(n);

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

	[~,VR,rep] = var2riss(A,V,y,r);
    if sserror(rep,y), continue; end % check DARE report, bail out on error

    DVR = diag(VR);
    F(r,y) = log(DVR)-LDV(r);

    if stats1
        stats.single.LR.tstat(r,y) = F(r,y);            % likelihood-ratio test statistic
        stats.single.F.tstat(r,y)  = DVR./DV(r) - 1;    % F-test statistic
    end

    if stats2
		[~,VR] = tsdata_to_var(X(r,:,:),p,regmode);
		DVR = diag(VR);
        stats.dual.LR.tstat(r,y) = log(DVR) - LDV(r); % likelihood-ratio test statistic
        stats.dual.F.tstat(r,y)  = DVR./DV(r) - 1;    % F-test statistic
    end
end

if stats1
    stats.single.F.pval  = 1-fcdf(K*stats.single.F.tstat,d,d2);
    stats.single.LR.pval = 1-chi2cdf(M*stats.single.LR.tstat,d);
end

if stats2
    stats.dual.F.pval    = 1-fcdf(K*stats.dual.F.tstat,d,d2);
    stats.dual.LR.pval   = 1-chi2cdf(M*stats.dual.LR.tstat,d);
end
