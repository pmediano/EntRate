% Convert covariance matrix to correlation matrix
%
% Result guaranteed positive-definite if covariance matrix is
% positive-definite (it should be! - check p == 0), otherwise
% a result is still returned, but should be treated as suspect.

function [R,p] = cov2corr(V)

assert(ismatrix(V),             'covariance matrix must be a square matrix');
assert(size(V,1) == (size(V,2)),'covariance matrix must be square');

[R,p] = chol(V); % right Cholesky factor
if p > 0 % fall back on "safe" method (doesn't ensure positive-definite result)
	d = 1./sqrt(diag(V));
	R = bsxfun(@times,bsxfun(@times,d,V),d')
else      % ensures positive-definite result
	R = bsxfun(@rdivide,R,sqrt(sum(R.*R)));
	R = R'*R;
end
