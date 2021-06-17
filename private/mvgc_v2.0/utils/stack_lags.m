function [X0,XX] = stack_lags(X,p)

[n,m,N] = size(X);
assert(p < m,'too many lags');
p1 = p+1;

M = N*(m-p);

X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
XX = zeros(n,p,M);
for k = 1:p
    XX(:,k,:) = reshape(X(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
end
XX = reshape(XX,n*p,M);         % stack lags

if nargout < 2
    X0 = [X0;XX];
end
