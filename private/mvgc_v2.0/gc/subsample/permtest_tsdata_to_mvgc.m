%% permtest_tsdata_to_mvgc
%
% Calculate null distribution for conditional time-domain MVGC from time series
% data, based on a permutation test
%
% <matlab:open('permtest_tsdata_to_mvgc.m') code>
%
%% Syntax
%
%     FP = permtest_tsdata_to_mvgc(U,x,y,p,bsize,nsamps,regmode,acmaxlags,acdectol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     p          model order (number of lags)
%     bsize      permutation block size (default: use model order)
%     nsamps     number of permutations
%     regmode    regression mode (default as for 'tsdata_to_var')
%     acmaxlags  maximum autocovariance lags  (default as for 'var_to_autocov')
%     acdectol   autocovariance decay tolerance (default as for 'var_to_autocov')
%
% _output_
%
%     FP         permutation test Granger causalities (null distribution)
%
%% Description
%
% Returns |nsamps| samples from the empirical null distribution of the
% time-domain MVGC from the variable |Y| (specified by the vector of indices
% |y|) to the variable |X| (specified by the vector of indices |x|), conditional
% on all other variables in the time series data |U|, based on randomly
% permuting blocks of size |bsize| of the source variable |Y| [2]. |p| is the
% model order; for other parameters see <tsdata_to_var.html |tsdata_to_var|> and
% <var_to_autocov.html |var_to_autocov|>.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] M. J. Anderson and J. Robinson, Permutation tests for linear models,
% _Aust. N. Z. J. Stat._ 43(1), 2001.
%
%% See also
%
% <mvgc_demo_permtest.html |mvgc_demo_permtest|> |
% <permtest_tsdata_to_pwcgc.html |permtest_tsdata_to_pwcgc|> |
% <permtest_tsdata_to_smvgc.html |permtest_tsdata_to_smvgc|> |
% <permtest_tsdata_to_spwcgc.html |permtest_tsdata_to_spwcgc|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_mvgc.html |autocov_to_mvgc|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function FP = permtest_tsdata_to_mvgc(U,x,y,p,bsize,nsamps,regmode,~,~) % v2.0 - last two parameters no longer required

if nargin < 7, regmode   = []; end % ensure default

if isempty(bsize), bsize = p; end % default to model order

[n,m,N] = size(U);
assert(m > p,'too many lags');

x = x(:)'; % vectorise
y = y(:)'; % vectorise
assert(all(x >=1 & x <= n),     'some x indices out of range');
assert(all(y >=1 & y <= n),     'some y indices out of range');
assert(isempty(intersect(x,y)), 'x and y indices must be distinct');
ny = length(y);

assert(isscalar(bsize) && isint(bsize) && bsize > 0, 'block size must be a positive integer');
nblocks  = floor(m/bsize); % number of blocks

if nblocks*bsize ~= m
    oldm = m;
    m = nblocks*bsize;
    U = U(:,1:m,:);
    fprintf(2,'WARNING: truncating sequence length by %d observations\n',oldm-m);
end

FP = nan(nsamps,1);

Y = reshape(U(y,:,:),ny,bsize,nblocks,N); % stack blocks of y-variable

for s = 1:nsamps
    fprintf('MVGC: permutation test sample %d of %d',s,nsamps);
    
    % generate permutation time series: "block permute" the y-variable per-trial
    
    for r = 1:N
        U(y,:,r) = reshape(Y(:,:,randperm(nblocks),r),ny,m); % permute blocks and unstack
    end
    
    % estimate permutation test VAR parameters
    
    [AP,SIGP] = tsdata_to_var(U,p,regmode);
    if isbad(AP), fprintf(2,' *** VAR estimation failed\n'); continue; end % something went badly wrong
    
    % calculate permutation test MVGC
    
    [AP1,CP,KP,res] = var_to_ss(AP,SIGP,2);
    if res.error, fprintf(2,' *** %s\n',res.errmsg); continue; end
    FP(s) = ss_to_mvgc(AP1,CP,KP,SIGP,x,y);
    
    fprintf('\n');
end
