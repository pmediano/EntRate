% Lag-stacking routine that can handle multiple trials of varying length

function XX = stacklags(X,p)

p1 = p+1;

if iscell(X)

	assert(isvector(X),'X must be a matrix, a 3-dim array, or a cell vector of matrices');
	n = size(X{1},1);
	N = length(X);
	m = zeros(N,1);
	for r = 1:N
		assert(ismatrix(X{r}),'trials data must be matrices');
		[n1,m(r)] = size(X{r});
		assert(n1 == n,'number of variables in trials don''t match');
	end
	assert(all(p < m),'too many lags');

	o = [0;cumsum(m-p)]; % trial length offsets
	M = o(end);

	XX = zeros(n,p1,M);
	for k = 0:p % concatenate trials for k-lagged observations
		for r = 1:N
			XX(:,k+1,o(r)+1:o(r+1)) = X{r}(:,p1-k:m(r)-k);
		end
	end

else

	assert(ismatrix(X) || ndims(X) == 3,'X must be a matrix, a 3-dim array, or a cell vector of matrices');
	[n,m,N] = size(X);
	assert(p < m,'too many lags');

	M = N*(m-p);

	XX = zeros(n,p1,M);
	for k = 0:p % concatenate trials for k-lagged observations
		XX(:,k+1,:) = reshape(X(:,p1-k:m-k,:),n,M);
	end

end

XX = reshape(XX,n*p1,M); % stack lags
