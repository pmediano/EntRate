% Version of 'demean' that can handle multiple trials of varying length

function [X,xmean,xstd] = demean_mt(X,normalise)

if nargin < 2 || isempty(normalise), normalise = false; end

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

	o = [0;cumsum(m)]; % trial length offsets
	M = o(end);

	Y = zeros(n,M);
	for r = 1:N
		Y(:,o(r)+1:o(r+1)) = X{r};
	end
	xmean = mean(Y,2);
	Y  = bsxfun(@minus,Y,xmean);
	if normalise
		xstd = std(Y,[],2);
		Y = bsxfun(@rdivide,Y,xstd);
	elseif nargout > 2
		xstd = std(Y,[],2);
	end
	for r = 1:N
		X{r} = Y(:,o(r)+1:o(r+1));
	end

else

	assert(ismatrix(X) || ndims(X) == 3,'X must be a matrix, a 3-dim array, or a cell vector of matrices');
	[n,m,N] = size(X);
	X = X(:,:);
	xmean = mean(X,2);
	X = bsxfun(@minus,X,xmean);
	if normalise
		xstd = std(X,[],2);
		X = bsxfun(@rdivide,X,xstd);
	elseif nargout > 2
		xstd = std(X,[],2);
	end
	X = reshape(X,n,m,N);

end
