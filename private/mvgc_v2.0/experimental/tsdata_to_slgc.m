function [FLR,stats] = tsdata_to_slgc(U,x,y,p,verb)

% Single-lag GC of y -> x at 1:p lags

if nargin < 5 || isempty(verb), verb = false; end

[n,m,N] = size(U);
assert(N == 1,'sorry, multi-trial not yet implemented')
assert(p < m,'too many lags');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];               % indices of reduced variables

nx = length(x);
ny = length(y);
nr = n-ny;

M  = N*(m-p); % chi^2 scaling factor = effective number of observations
d  = nx*ny;   % chi^2 df and F df1
stats.LR.mean  = d/M;
stats.LR.cdf   = @(t) chi2cdf(M*t,d);
stats.LR.icdf  = @(p) chi2inv(p,d)/M;

d2 = nx*(M-n)-1; % F df2
K  = d2/d;       % F scaling factor
stats.F.mean  = d/(d2-2);
stats.F.cdf   = @(t) fcdf(K*t,d,d2);
stats.F.icdf  = @(p) finv(p,d,d2)/K;

p1   = p+1;
pn   = p*n;
pnr  = p*nr;
pny1 = (p-1)*ny;

X0 = U(x,p1:m);
UL = zeros(n,p,M);
for k = 1:p
	UL(:,k,:) = reshape(U(:,p1-k:m-k,:),n,M);
end
UL = reshape(UL,pn,M);
E = X0-(X0/UL)*UL;
Vx = (E*E')/(M-1);
LDVx = logdet(Vx);

FLR = nan(p,1);
FF  = nan(p,1);

UL = reshape(UL,[n p M]);
XZL = reshape(UL(r,:,:),pnr,M);
for k = 1:p
	if verb, fprintf('lag = %2d of %2d\n',k,p); end
	ULR = [XZL;reshape(UL(y,[1:k-1 k+1:p],:),pny1,M)]; % omit k-th y lag
	ER = X0-(X0/ULR)*ULR;
	VRx = (ER*ER')/(M-1);
	FLR(k) = logdet(VRx) - LDVx;        % log-likelihood-ratio test statistic
	FF(k)  = trace(VRx)/trace(Vx) - 1;  % F-test statistic
end

stats.LR.tstat = FLR;
stats.F.tstat  = FF;
