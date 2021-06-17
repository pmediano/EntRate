function [LL,Lk,Lm] = tsdata_to_mle_mt(X,q,verb)

vlt = iscell(X);

demean_mt(X,false);

q1 = q+1;

if vlt
	assert(isvector(X),'X must be a matrix, a 3-dim array, or a cell vector of matrices');
	n = size(X{1},1);
	N = length(X);
	m = zeros(N,1);
	for r = 1:N
		assert(ismatrix(X{r}),'trials data must be matrices');
		[n1,m(r)] = size(X{r});
		assert(n1 == n,'number of variables in trials don''t match');
	end
	assert(all(q < m),'too many lags');
else
	assert(ismatrix(X) || ndims(X) == 3,'X must be a matrix, a 3-dim array, or a cell vector of matrices');
	[n,m,N] = size(X);
	assert(q < m,'too many lags');
end

% Note: p = 1 is order 0 term!
LL = nan(q1,1); % log-likelihood
Lk = nan(q1,1); % number free parameters
Lm = nan(q1,1); % effective sample size

% order zero likelihood

if vlt
	o = [0;cumsum(m)];          % trial length - model order offsets
	M = o(end);
	X0 = zeros(n,M);
	for r = 1:N                 % concatenate trials for unlagged observations
		X0(:,o(r)+1:o(r+1)) = X{r}(:,q1:m(r));
	end
else
	M  = N*m;
	X0 = reshape(X,n,M);        % concatenate trials for unlagged observations
end
DSIG = det((X0*X0')/(M-1)); % covariance matrix determinant
assert(DSIG > 0,'covariance matrix not positive-definite');
LL(1) = -(M/2)*log(DSIG);
Lk(1) =  0;
Lm(1) =  M;

% loop through model orders

for p = 1:q

	if verb, fprintf('model order = %d',p); end

	p1 = p+1;

	if vlt
		o = [0;cumsum(m-p)];           % trial length - model order offsets
		M = o(end);
		X0 = zeros(n,M);
		for r = 1:N                    % concatenate trials for unlagged observations
			X0(:,o(r)+1:o(r+1)) = X{r}(:,p1:m(r));
		end
		XL = zeros(n,p,M);
		for k = 1:p                    % concatenate trials for k-lagged observations
			for r = 1:N
				XL(:,k,o(r)+1:o(r+1)) = X{r}(:,p1-k:m(r)-k);
			end
		end
	else
		M  = N*(m-p);
		X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
		XL = zeros(n,p,M);
		for r = 1:p
			XL(:,r,:) = reshape(X(:,p1-r:m-r,:),n,M); % concatenate trials for p-lagged observations
		end
	end
	XL = reshape(XL,n*p,M);         % stack lags

	wstate = warning('off','all'); lastwarn('');
	A = X0/XL;                     % OLS using QR decomposition
	wmsg = lastwarn; warning(wstate);
	if ~isempty(wmsg) % rank-deficient?
		if ~verb, fprintf('model order = %d',p); end
		fprintf(2,'  WARNING: VAR estimation may be problematic (%s)',wmsg);
		% not necessarily a show-stopper - carry on
	end
	if isbad(A)                     % something went badly wrong
		if ~verb, fprintf('model order = %d',p); end
		fprintf(2,'  WARNING: VAR estimation failed\n');
		continue % show-stopper
	end

	E    = X0-A*XL;                % residuals
	DSIG = det((E*E')/(M-1));      % residuals covariance matrix determinant

	if DSIG <= 0
		if ~verb, fprintf('model order = %d',p); end
		fprintf(2,'  WARNING: residuals covariance not positive definite\n');
		continue % show-stopper
	end

	LL(p+1) = -(M/2)*log(DSIG);
	Lk(p+1) =  p*n*n;
	Lm(p+1) =  M;

	if verb, fprintf(1,'\n'); end
end
