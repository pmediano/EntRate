function [S,K,J] = mardia(X,debias)

% Mardia's multivariate skew/kurtosis and Koizumi-Okamoto-Seo's multivariate
% Jarque-Bera test statistics.

if nargin < 2 || isempty(debias), debias = false; end

S = NaN;
K = NaN;

X = X(:,:);
[n,m] = size(X);

d = n*(n+2);

X  = bsxfun(@minus,X,mean(X,2)); % demean

V = (X*X')/m;

[U,cholp] = chol(V,'lower');
if cholp, return; end % caller must test for NaN

w = U\X;
W = w'*w;

W3 = W.^3;
S = mean(W3(:));          % Mardia's multivariate skewness: for n = 1, S = skewness(n)^2

K = mean(diag(W).^2) - d; % Mardia's multivariate kurtosis: for n = 1, K = kurtosis(n)-3

if debias
    S = S - d*((m+1)*(n+1)-6)/((m+1)*(m+3));
    K = K + 2*d/(m+1);
end

if nargout > 2
    J = S/6 + (K^2)/(8*d); % multivariate Jarque-Bera test statistic [Koizumi, Okamoto and  Seo, J. Stat. Adv. Theo. Appl. 1(2), 2009]
end
