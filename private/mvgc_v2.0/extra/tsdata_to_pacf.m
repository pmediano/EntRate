function [R,M] = tsdata_to_pacf(X,p,fast)

% Compute partial autocorrelation function to p lags for variables in X

if nargin < 3 || isempty(fast), fast = true; end

[n,m,N] = size(X);
assert(p < m,'too many lags');

assert(isscalar(fast) && (islogical(fast) || isint(fast)),'Flag must convert to a logical scalar');

X = demean(X);

R = nan(n,n,p);

wstate = warning('off','all');
if fast

    for k = 1:p

        M = N*(m-k);
        X0 = reshape(X(:,k+1:m,:),n,M);
        XL = zeros(n,k,M);
        for r = 1:k
            XL(:,r,:) = reshape(X(:,k+1-r:m-r,:),n,M);
        end
        XL = reshape(XL,n*k,M);

        A = X0/XL;
        A = reshape(A,n,n,k);

        R(:,:,k) = A(:,:,k);

    end

else

    for k = 1:p

        M = N*(m-k);
        X0 = reshape(X(:,k+1:m,:),n,M);
        XL = zeros(n,k,M);
        for r = 1:k
            XL(:,r,:) = reshape(X(:,k+1-r:m-r,:),n,M);
        end

        Xp = reshape(XL(:,k,:),n,M);
        XP = reshape(XL(:,1:k-1,:),n*(k-1),M);

        E0 = X0-(X0/XP)*XP;
        Ep = Xp-(Xp/XP)*XP;

        R(:,:,k) = corr(E0',Ep');

    end

end
warning(wstate);

M = N*(m-(1:p)');

%alpha = alpha/(n*n*p);              % Bonferroni correction:
%Rcrit = norminv(1-alpha/2)/sqrt(M); % 2-sided
