function sig = ssignificance(pval,alpha,correction,sym)

if ~sym
    sig = significance(pval,alpha,correction);
    return
end

% Like 'significance()', but assumes a symmetric matrix and uses only upper triangle (diagonals are ignored).

[n,n1] = size(pval);
assert(ismatrix(pval) && n1 == n,'p-values must be a square matrix');

utidx = logical(triu(ones(n),1)); % logical indices of upper triangle

utsig = significance(pval(utidx),alpha,correction);

sig = NaN(n);
sig(utidx) = utsig;
sig = symmetrise(sig);
