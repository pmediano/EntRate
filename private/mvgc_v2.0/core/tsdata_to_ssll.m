function [LL,Lk,Lm,sval,SVC,r] = tsdata_to_ssll(X,h,r,verb)

% NOTE: returns *AVERAGE* log-likelihood

[n,m,N] = size(X);

assert(all(isint(h(:))),'past/future horizon must be a 2-vector or a scalar positive integer');
if isscalar(h)
    p = h;    f = h;
elseif isvector(h) && length(h) == 2
    p = h(1); f = h(2);
else
    error('past/future horizon must be a 2-vector or a scalar positive integer');
end
assert(p+f < m,'past/future horizon too large (or not enough data)');
rmax = n*min(p,f);

if nargin < 3 || isempty(r),    r      = rmax;  end
if nargin < 4 || isempty(verb), verb   = false; end

assert(isscalar(r) && isint(r) && r >= 0 && r <= rmax,'model order must be empty, or a positive integer <= n*min(p,f) = %d',rmax);

X = demean(X); % no constant term (don't normalise!)

LL  = nan(r,1); % log-likelihood
Lk  = nan(r,1); % number free parameters
Lm  = nan(r,1); % effective sample size
SVC = nan(r,1); % effective sample size

mp  = m-p;
mp1 = mp+1;
mf  = m-f;
mh  = mp1-f; % m-p-f+1

M  = N*mp;
M1 = N*mp1;
Mh = N*mh;

Xf = zeros(n,f,mh,N);
for k = 1:f
    Xf(:,k,:,:) = X(:,p+k:mf+k,:);
end
Xf = reshape(Xf,n*f,Mh);

XP = zeros(n,p,mp1,N);
for k = 0:p-1
    XP(:,k+1,:,:) = X(:,p-k:m-k,:);
end
Xp = reshape(XP(:,:,1:mh,:),n*p,Mh);
XP = reshape(XP,n*p,M1);

[Wf,cholp] = chol((Xf*Xf')/Mh,'lower');
assert(cholp == 0,'forward weight matrix not positive definite');

[Wp,cholp] = chol((Xp*Xp')/Mh,'lower');
assert(cholp == 0,'backward weight matrix not positive definite');

BETA = Xf/Xp; % 'OH' estimate: regress future on past
assert(all(isfinite(BETA(:))),'subspace regression failed');

[~,S,U] = svd(Wf\BETA*Wp); % SVD of CCA-weighted OH estimate

sval = diag(S);    % the singular values
Lk   = 2*n*(1:r)'; % number of free parameters (Hannan & Deistler, see also Bauer 2001) ... or r*r+2*n*r ???
SVC  = -log(1-[sval(2:end);0]) + Lk*(log(Mh)/Mh); % Bauer's Singular Value Criterion
Lm   =  M*ones(r,1);

XX = reshape(X(:,p+1:m,:),n,M);
ssval = sqrt(sval);

for k = 1:r

	if verb, fprintf('model order = %d',k); end

	Z = reshape((diag(ssval(1:k))*U(:,1:k)'/Wp)*XP,k,mp1,N); % Kalman states estimate; note that Z starts at t = p+1, has length mp1 = m-p+1
	ZZ = reshape(Z(:,1:mp,:),k,M);

	C = XX/ZZ;        % observation matrix
	if ~all(isfinite(C(:)))
		if ~verb, fprintf('model order = %d',k); end
		fprintf(2,'  WARNING: observation matrix estimation failed\n');
		continue;
	end

	E = XX - C*ZZ; % innovations
	V = (E*E')/M;  % innovations covariance matrix (ML estimator)

	LL(k) = -logdet(V)/2; % approximation to avergae log-likelihood

	if verb, fprintf('\n'); end

end
