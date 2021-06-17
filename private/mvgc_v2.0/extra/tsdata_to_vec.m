function [r,h,pval,stat,cval,eval,evec,rLL,uLL] = tsdata_to_vec(y,cmode,p,alpha,summary)

if nargin < 5 || isempty(summary), summary = true; end
if summary, distr = 'summary'; else distr = 'off'; end

% IMPORTANT NOTE: The 'lags' parameter to 'jcitest' is the VEC autoregressive
% order = (VAR order)-1 !!!

warning('off','econ:jcitest:RightTailStatTooBig');
warning('off','econ:jcitest:RightTailStatTooSmall');

[h,pval,stat,cval,mles] = jcitest(y','model',cmode,'lags',p,'test',{'trace','maxeig'},'alpha',alpha,'display',distr);
warning('on','econ:jcitest:RightTailStatTooBig');
warning('on','econ:jcitest:RightTailStatTooSmall');

% get rid of those silly tables

h    = h{:,:}';    % 'true' indicates rejection of the null of cointegration rank r in favor of the alternative
pval = pval{:,:}';
stat = stat{:,:}';
cval = cval{:,:}';

n = size(y,1);
r = zeros(1,2); % r = 0 => unit root but not cointegrated, r = n => stable VAR
for i = 1:2 % trace, maxeig
    rr = find(~h(:,i),1)-1; % working from the left, number of 'true's before we hit a 'false'
    if isempty(rr), r(i) = n; else r(i) = rr; end
end

if nargout > 5
    eval = zeros(n,1);
    evec = zeros(n,n);
    for k = 1:n
        eval(k)   = mles{1,k}.eigVal;
        evec(:,k) = mles{1,k}.eigVec;
    end
end

if nargout > 7
    rLL = zeros(n,1);
    for k = 1:n
        rLL(k) = mles{1,k}.rLL;
    end
    uLL = zeros(n,2);
    for i = 1:2 % trace, maxeig
        for k = 1:n
            uLL(k,i) = mles{i,k}.uLL;
        end
    end
end
