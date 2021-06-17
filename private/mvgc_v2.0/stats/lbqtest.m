function Q = lbqtest(X,p,hmax,standardised)

% Ljung-Box Q "portmanteau" test statistic for serial correlation.
%
% Rule-of-thumb: set hmax ~ log(m). If testing residuals of an AR model, you
% might set hmax to the largest model order tested for model selection (e.g. by
% AIC, BIC, etc.), which you (should) choose on domain-specific grounds.

if nargin < 4 || isempty(standardised), standardised = false; end

assert(hmax > p,'Maximum autocorrelation lag must be larger than model order');

[n,m,N] = size(X);

assert(hmax < m,'Maximum autocorrelation lag must be smaller than number of observations');

T = N*m;

C = zeros(n);
for i = 1:N
    C = C + (X(:,:,i)*X(:,:,i)')/m;
end
C = C/N; % covariance matrix

I = C\eye(n); % inverse covariance matrix

Q = zeros(hmax,1);
for h = 1:hmax
    C = zeros(n);
    for i = 1:N
        C = C + (X(:,h+1:m,i)*X(:,1:m-h,i)')/m;
    end
    C = C/N; % autcovariance matrix at lag h
    Q(h) = trace(C*I*C'*I)/(T-h);
end
Q = T*(T+2)*cumsum(Q);

Q(1:p) = NaN; % undefined for h <= p

if standardised
    df = n*n*(1:hmax-p)'; % chi2 degrees of freedom
    Q(p+1:h) = (Q(p+1:h)-df)./sqrt(2*df);
end
