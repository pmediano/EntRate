function [F,stats] = var_to_gwiogc(A,V,group,inout,whichstats,X,regmode)

% in/out group GCs

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(group,n);

switch lower(inout)
	case 'in',  gcin = true;
	case 'out', gcin = false;
	otherwise, error('in/out parameter must be ''in'' or ''out''');
end

stats1 = false;
stats2 = false;
if nargin > 4 && ~isempty(whichstats)
	switch lower(whichstats)
		case 'single', stats1 = true;
		case 'dual',   stats2 = true;
		case 'both',   stats1 = true; stats2 = true;
		case 'none'    % do nothing
		otherwise error('bad stats specification');
	end
end

stats = [];
if stats1 || stats2;
	if     stats2
		assert(nargin > 5 && ~isempty(X) && ~isscalar(X),'must supply time series data for dual-regression stats');
		assert(nargin > 6,'must supply regression mode (same as for parameter estimates) for dual-regression stats');
        [~,m,N] = size(X);
	elseif stats1
		assert(nargin > 5 && ~isempty(X),'must supply number of observations and trials for single-regression stats');
		if isscalar(X)
			assert(nargin > 4 && ~isempty(regmode),'must supply number of observations and trials for single-regression stats');
			m = X;
			N = regmode;
			clear X regmode;
		else
			[~,m,N] = size(X);
		end
	end
    M  = N*(m-p); % chi^2 scaling factor = effective number of observations
    d  = nan(g,1);  % chi^2 df and F df1
    d2 = nan(g,1);  % F df2
    for a = 1:g
		if gcin
			nx = length(group{a});
			ny = n-nx;
		else
			ny = length(group{a});
			nx = n-ny;
		end
		d(a)  = p*nx*ny;
		d2(a) = nx*(M-p*n)-1;
		K  = d2./d;
		Xm = d/M;
		Fm = d./(d2-2);
	end
    statsf = struct('tstat',nan(g,1),'pval',nan(g,1));
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

F = nan(g,1);

for a = 1:g
	if gcin
		x = group{a};
		y = 1:n; y(x) = [];
	else
		y = group{a};
		x = 1:n; x(y) = [];
	end

	[~,VR,rep] = var2riss(A,V,y,x);
	if sserror(rep), return; end % check DARE report, bail out on error

	F(a) = logdet(VR) - logdet(V(x,x));

	if stats1
		stats.single.F.tstat(a)  = trace(VR)/trace(V(x,x)) - 1; % F-test statistic
		stats.single.LR.tstat(a) = F(a);                        % likelihood-ratio test statistic
	end

	if stats2
		[~,VR2] = tsdata_to_var(X(x,:,:),p,regmode);
		stats.dual.F.tstat(a)  = trace(VR2)/trace(V(x,x)) - 1;  % F-test statistic
		stats.dual.LR.tstat(a) = logdet(VR2) - logdet(V(x,x));  % likelihood-ratio test statistic
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
