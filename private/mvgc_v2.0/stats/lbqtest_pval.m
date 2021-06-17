function pval = lbqtest_pval(Q,n,p,h)

% p-values values for Ljung-Box Q "portmanteau" test statistic for serial correlation.

assert(isequal(size(Q),size(h)),'Q statistics don''t match autocorrelation lengths');

h(h <= p) = NaN;
pval = 1-chi2cdf(Q,n*n*(h-p));
