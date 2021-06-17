function [Y,nouts] = routl_m(X,madfac,repmed)

% Replace outliers in each variable in matrix X (variables x observations x
% trials), by a random non-outlier, or median if repmed is set. Outliers are
% defined as less than (median - madfac*mad) or greater than (median +
% madfac*mad) where mad is the median absoulute deviation of X.

assert(madfac> 0);

madfac = norminv(3/4)*madfac;

[n,m,N] = size(X);
if N > 1 % multi-trial
    X = reshape(X,n,N*m);
end

Y = X;
nouts = zeros(n,1);

for i = 1:n
    medi     = median(X(i,:));
    madi     = madfac*mad(X(i,:),1);
    outs     = (X(i,:) < medi-madi) | (X(i,:) > medi+madi);
    nouts(i) = nnz(outs);
    Z        = X(i,~outs); % non-outliers
    if repmed
        Y(i,outs) = median(Z)*ones(1,nouts(i));       % replace with median of non-outliers
    else
        Y(i,outs) = Z(randi(length(Z),[1 nouts(i)])); % replace with random non-outlier
    end
end

if N > 1 % multi-trial
    Y = reshape(Y,[n m N]);
end
