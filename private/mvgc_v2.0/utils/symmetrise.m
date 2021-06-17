function y = symmetrise(x,hermitian,ut2lt)

% Symmetrise a square matrix by reflecting the upper (default) or lower triangle
% across the diagonal. The default is to make a Hermitian matrix, which requires
% that the diagonal is real (if the entire matrix is real then the 'hermitian'
% flag has no effect). The diagonal is always left unchanged.

if nargin < 2 || isempty(hermitian), hermitian = true; end
if nargin < 3 || isempty(ut2lt),     ut2lt     = true; end

[n,n1] = size(x);
assert(ismatrix(x) && n1 == n,'Input must be a square matrix');

if ut2lt
    tidx = logical(tril(ones(n),-1)); % logical indices of lower triangle
else
    tidx = logical(triu(ones(n),+1)); % logical indices of upper triangle
end

y = x;
if hermitian
    assert(isreal(diag(x)),'Diagonal must be real for a Hermitian matrix');
    x = x';
else
    x = x.';
end
y(tidx) = x(tidx);
