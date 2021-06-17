function [R,I] = pwparcorr(V)

% Calculate pairwise partial correlation and (optionally) conditional mutual information matrices

[n,n1] = size(V);
assert(ismatrix(V) && n1 == n,'covariance matrix must be square');

R = NaN(n);
I = NaN(n);

[V cholp] = chol(V);
if cholp ~= 0, return; end  % error: not pos-def - user must check
V = inv(V);
R = bsxfun(@times,1./sqrt(sum(V.*conj(V),2)),V);
R = -R*R';         % NOTE: This means that -R is positive-definite, but ...
R(1:n+1:n*n) = 1;  % put 1's on the diagonal, so that "each variable partially-correlates perfectly with itself"; now 2I-R is positive-definite
if nargout > 1
    I = -log(1-R.*conj(R)); % will have Inf's on diagonal
end
