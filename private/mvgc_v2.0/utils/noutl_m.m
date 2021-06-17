function nouts = routl_m(X,madfac)

% Count outliers. Outliers are defined as less than (median - madfac*mad) or
% greater than (median + madfac*mad) where mad is the median absoulute deviation
% of X.

if madfac < eps % do nothing
    nouts = 0;
    return;
end

madfac = norminv(3/4)*madfac;

[n,m,N] = size(X);
X = X(:,:);

nouts = zeros(n,1);
for i = 1:n
    medi     = median(X(i,:));
    madi     = madfac*mad(X(i,:),1);
    outs     = X(i,:) < medi-madi | X(i,:) > medi+madi;
    nouts(i) = nnz(outs);
end
