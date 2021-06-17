function [FLR,stats] = tsdata_to_pwslgc(U,p,verb)

% Single-lag GC of y -> x at 1:p lags

if nargin < 3 || isempty(verb), verb = false; end

[n,m,N] = size(U);
assert(N == 1,'sorry, multi-trial not yet implemented')
assert(p < m,'too many lags');

nr = n-1;

M  = N*(m-p); % chi^2 scaling factor = effective number of observations
d  = 1;       % chi^2 df and F df1
stats.LR.mean  = d/M;
stats.LR.cdf   = @(t) chi2cdf(M*t,d);
stats.LR.icdf  = @(p) chi2inv(p,d)/M;

d2 = M-n-1;   % F df2
K  = d2/d;    % F scaling factor
stats.F.mean  = d/(d2-2);
stats.F.cdf   = @(t) fcdf(K*t,d,d2);
stats.F.icdf  = @(p) finv(p,d,d2)/K;

p1   = p+1;
pn   = p*n;
pnr  = p*nr;
pny1 = p-1;

U0 = U(:,p1:m);
UL = zeros(n,p,M);
for k = 1:p
	UL(:,k,:) = reshape(U(:,p1-k:m-k,:),n,M);
end
UL = reshape(UL,pn,M);
E = U0-(U0/UL)*UL;
DV = sum(E.*E,2)/(M-1);
LDV = log(DV);

FLR = nan(p,n,n);
FF  = nan(p,n,n);

UL = reshape(UL,[n p M]);

for y = 1:n
	if verb, fprintf('source = %2d of %2d ',y,n); end
	r = [1:y-1,y+1:n];
	UR = reshape(UL(r,:,:),pnr,M);
	for k = 1:p
		if verb, fprintf('.'); end
		ULR = [UR;reshape(UL(y,[1:k-1 k+1:p],:),pny1,M)]; % omit k-th y lag
		ER = U0(r,:)-(U0(r,:)/ULR)*ULR;
		DVR = sum(ER.*ER,2)/(M-1);
		FLR(k,r,y) = log(DVR) - LDV(r); % likelihood-ratio test statistic
		FF (k,r,y) = DVR./DV(r) - 1;    % F-test statistic
	end
	if verb, fprintf('\n'); end
end

stats.LR.tstat = FLR;
stats.F.tstat  = FF;
