% Version of 'tsdata_to_var' that can handle multiple trials of varying length (OLS only)

function [A,SIG,E] = tsdata_to_var_mt(X,p,dm)

if nargin < 3 || isempty(dm), dm = true; end

if dm, demean_mt(X,false); end % no constant term!

vlt = iscell(X);

p1 = p+1;

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
	assert(all(p < m),'too many lags');

	% store lags

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

	assert(ismatrix(X) || ndims(X) == 3,'X must be a matrix, a 3-dim array, or a cell vector of matrices');
	[n,m,N] = size(X);
	assert(p < m,'too many lags');

	% store lags

    M = N*(m-p);
    X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
    XL = zeros(n,p,M);
    for k = 1:p                    % concatenate trials for k-lagged observations
        XL(:,k,:) = reshape(X(:,p1-k:m-k,:),n,M);
    end

end

XL = reshape(XL,n*p,M);            % stack lags

% do OLS

A = X0/XL;                         % OLS (via QR decomposition)
if ~all(isfinite(A(:)))            % something went badly wrong
	return
end

if nargout > 1

	% calculate residuals

	if vlt
        EE  = X0-A*XL;             % residuals
        SIG = (EE*EE')/(M-1);      % residuals covariance matrix
        E   = cell(size(X));       % put residuals back into per-trial form
		for r = 1:N
			E{r} = EE(:,o(r)+1:o(r+1));
		end
	else
        E   = X0-A*XL;             % residuals
        SIG = (E*E')/(M-1);        % residuals covariance matrix
        E   = reshape(E,n,m-p,N);  % put residuals back into per-trial form
	end
end

A = reshape(A,n,n,p);              % so A(:,:,k) is the k-lag coefficients matrix
