function [F,stats] = var_to_mvgc(A,V,x,y,whichstats,X,regmode)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];               % indices of reduced variables

nx = length(x);
ny = length(y);
xr = 1:nx;               % index of x in reduced variables

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
if stats1 | stats2;
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
    M  = N*(m-p);      % chi^2 scaling factor = effective number of observations
    d  = p*nx*ny;      % chi^2 df and F df1
    Xm = d/M;          % chi^2 mean
    d2 = nx*(M-p*n)-1; % F df2
    K  = d2/d;         % F scaling factor
    Fm = d/(d2-2);     % F mean
    statsf = struct('tstat',nan,'pval',nan);
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

F = NaN;

[~,VR,rep] = var2riss(A,V,y,r);
if sserror(rep), return; end % check DARE report, bail out on error

F = logdet(VR(xr,xr)) - logdet(V(x,x));

if stats1
	stats.single.F.tstat  = trace(VR(xr,xr))/trace(V(x,x)) - 1;  % F-test statistic
	stats.single.F.pval   = 1-fcdf(K*stats.single.F.tstat,d,d2);
	stats.single.LR.tstat = F;                                   % likelihood-ratio test statistic
	stats.single.LR.pval  = 1-chi2cdf(M*stats.single.LR.tstat,d);
end

if stats2
	[~,VR] = tsdata_to_var(X(r,:,:),p,regmode);
	stats.dual.F.tstat  = trace(VR(xr,xr))/trace(V(x,x)) - 1;    % F-test statistic
	stats.dual.F.pval   = 1-fcdf(K*stats.dual.F.tstat,d,d2);
	stats.dual.LR.tstat = logdet(VR(xr,xr)) - logdet(V(x,x));    % likelihood-ratio test statistic
	stats.dual.LR.pval  = 1-chi2cdf(M*stats.dual.LR.tstat,d);
end
