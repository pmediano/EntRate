function [F,stats] = var_to_gwggc(A,V,group,whichstats,X,regmode)

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
    d  = nan(g,1);  % chi^2 df and F df1
    d2 = nan(g,1);  % F df2
	for a = 1:g
		x = group{a};
		nx = length(x);
		ny = nx-1;
		d(a)  = p*nx*ny;
		d2(a) = nx*(M-p*n)-1;
	end
	K  = d2./d;
	Xm = d/M;
	Fm = d./(d2-2);
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

DV = diag(V);
LDV = log(DV);

for a = 1:g
	x = group{a};
    z = 1:n; z(x) = []; % omit x
    nx = length(x);

	Fa = 0;
	Ta  = 0;
	if stats1
		TRa = 0;
	end
	if stats2
		F2a  = 0;
		TR2a = 0;
	end
	sserr = false;
	for i = 1:nx
		xi = x(i);
		y = x; y(i) = [];
		r = [xi z];

		[~,VR,rep] = var2riss(A,V,y,r);
		if sserror(rep,i) % check DARE report, bail out on error
			sserr = true;
			break;
		end
		VRxi = VR(1,1);
		Fa = Fa + log(VRxi) - LDV(xi);
		Ta = Ta  + DV(xi);

		if stats1
			TRa = TRa + VRxi;
		end

		if stats2
			[~,VR2] = tsdata_to_var(X(r,:,:),p,regmode);
			VR2xi = VR2(1,1);
			F2a  = F2a  + log(VR2xi) - LDV(xi);
			TR2a = TR2a + VR2xi;
		end

	end
	if sserr, continue; end % bail out on earlier error

	F(a) = Fa;

	if stats1
		stats.single.F.tstat(a)  = TRa/Ta - 1;  % F-test statistic
		stats.single.LR.tstat(a) = Fa;          % likelihood-ratio test statistic
	end

	if stats2
		stats.dual.F.tstat(a)    = TR2a/Ta - 1; % F-test statistic
		stats.dual.LR.tstat(a)   = F2a;         % likelihood-ratio test statistic
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
