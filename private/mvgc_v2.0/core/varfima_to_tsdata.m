function [Y,X,mtrunc] = varfima_to_tsdata(A,B,fiparms,X,m,N,mtrunc,decfac)

% [Y,X,mtrunc] = varfima_to_tsdata(A,B,fiparms,V,m,N,mtrunc,decfac)
%
% or
%
% [Y,X,mtrunc] = varfima_to_tsdata(A,B,fiparms,X,[],[],mtrunc,decfac)

fint = ~isempty(fiparms);

if nargin < 7 || isempty(mtrunc) % automatic calculation - transients decay with rate given by VAR spectral radius
	assert(~fint,'Automatic truncation won''t work for fractional integration!');
	if isempty(B)
		mtb = 0;
	else
		if isvector(B), mtb = length(B); else, mtb = size(B,3); end
	end
	if isempty(A)
		mta = 0;
	else
		rho = specnorm(A);
		assert(rho < 1-eps,'Unstable: can''t auto-truncate');
		if nargin < 8 || isempty(decfac), decfac = 1; end
		mta = ceil(decfac*(-log(eps))/(-log(rho)));
	end
	mtrunc = mtb+mta;
else
    assert(isscalar(mtrunc) && isint(mtrunc) && mtrunc >= 0,'truncation parameter must be a non-negative integer');
end

[n,m1,N1] = size(X);
if nargin < 5 || isempty(m) % X is input
	assert(nargin < 6 || isempty(N),'number of trials specified by input!');
	m = m1;
	N = N1;
else          % X is a covariance matrix - generate Gaussian residuals as input
	if nargin < 6 || isempty(N), N = 1; end % single trial
	assert(ismatrix(X) && m1 == n,'covariance matrix not square');
	[L,cholp] = chol(X,'lower');
	assert(cholp == 0,'covariance matrix not positive-definite');
	m = m+mtrunc;
	X = zeros(n,m,N);
	for k = 1:N
		X(:,:,k) = L*randn(n,m); % Gaussian residuals
	end
end

assert(mtrunc < m,'too much truncation!');

if N > 1 % multi-trial
    Y = zeros(n,m,N);
	if fint
		for k = 1:N
			Y(:,:,k) = gdfilter(mvfilter(B,A,X(:,:,k)),fiparms);
		end
	else
		for k = 1:N
			Y(:,:,k) = mvfilter(B,A,X(:,:,k));
		end
	end
else
	if fint
		Y = gdfilter(mvfilter(B,A,X),fiparms);
	else
		Y = mvfilter(B,A,X);
	end
end

if mtrunc > 0
	Y = Y(:,mtrunc+1:m,:);
	if nargout > 1
		X = X(:,mtrunc+1:m,:);
	end
end
