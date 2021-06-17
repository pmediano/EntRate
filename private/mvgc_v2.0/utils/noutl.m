function nouts = noutl(X,sdfac)

% Count outliers. Outliers are defined as less than (mean - sdfac*sd) or greater
% than (mean + sdfac*sd) where sd is the standard deviation of X.

if sdfac < eps % do nothing
    nouts = 0;
    return;
end

[n,m,N] = size(X);
X = X(:,:);

nouts = zeros(n,1);
for i = 1:n
    xmean      = mean(X(i,:));
    xlim       = sdfac*std(X(i,:));
    outs       = X(i,:) < xmean-xlim | X(i,:) > xmean+xlim; % logical idx of outliers
    nouts(i)   = nnz(outs);  % number of outliers
end
