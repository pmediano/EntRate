function Q = lbqtest_cval(n,p,h,alpha)

% Critical values for Ljung-Box Q "portmanteau" test statistic for serial correlation.

h(h <= p) = NaN;
Q = chi2inv(1-alpha,n*n*(h-p));
