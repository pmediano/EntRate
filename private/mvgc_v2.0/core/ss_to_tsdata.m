function [X,Z,E,mtrunc] = ss_to_tsdata(A,C,K,V,m,N,mtrunc,decfac)

if nargin < 6 || isempty(N), N = 1; end % single trial

if nargin < 7 || isempty(mtrunc) % automatic calclation - transients decay with rate given by VAR spectral radius
	rhoa = max(abs(eig(A)));
	assert(rhoa < 1-eps,'Unstable: can''t auto-truncate');
	if nargin < 8 || isempty(decfac), decfac = 1; end
    mtrunc = ceil(decfac*(-log(eps))/(-log(rhoa)));
else
    assert(isscalar(mtrunc) && isint(mtrunc) && mtrunc >= 0,'truncation parameter must be a non-negative integer');
end

[n,r] = ss_parms(A,C,K,V);

% Generate Gaussian innovations with covariance matrix V

[VL,cholp] = chol(V,'lower');
assert(cholp == 0,'covariance matrix not positive-definite');

mtot = m+mtrunc;
E = nan(n,mtot+1,N);
Z = nan(r,mtot  ,N);
X = nan(n,mtot  ,N);
for u = 1:N
    E(:,:,u) = VL*randn(n,mtot+1);
    Z(:,:,u) = mvfilter([],A,K*E(:,1:mtot,u)); % state variable is AR(1) with innovations K*E (lagged)
    X(:,:,u) = C*Z(:,:,u) + E(:,2:mtot+1,u);   % observations variable uses unlagged innovations
end

E = E(:,2:mtot+1,:); % align innovations with X

if mtrunc > 0
    X = X(:,mtrunc+1:end,:);
    if nargout > 1
		Z = Z(:,mtrunc+1:end,:);
		if nargout > 2
			E = E(:,mtrunc+1:end,:);
		end
	end
end
