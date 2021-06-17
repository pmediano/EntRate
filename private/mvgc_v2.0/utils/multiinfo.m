function g = multiinfo(V,umean,iscorr)

% Return the multi-information of a (positive-definite) variance-covariance
% matrix V. This is equal to -log|R|, where R is the corresponding correlation
% matrix.
%
% If the flag 'umean' is set, then V must be an integer and we return -log<|R|>,
% where R is uniform random on the set of all correlation matrices of size V x V.
%
% If the 'iscorr' flag is set, return the correlation-like rho = sqrt(1-exp(-(2/n)*g))

if nargin < 2 || isempty(umean),  umean  = false; end
if nargin < 3 || isempty(iscorr), iscorr = false; end

if umean % return negative log of expected value of multi-information for uniform random correlation matrix
    assert(isscalar(V) && isint(V) && V > 1,'dimension must be an integer greater than 1');
    n = V;
	k = 1:n-1;
	g = sum(k.*(log(k+2)-log(k+1)));
else
	[n,n1] = size(V);
	assert(ismatrix(V) && n1 == n,'covariance matrix must be square');
	g = sum(log(diag(V)))-logdet(V);
end

if iscorr
	g = sqrt(1-exp(-(2/n)*g));
end
