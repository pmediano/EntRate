%% permtest_tsdata_to_pwcgc
%
% Calculate null distribution for pairwise-conditional time-domain MVGCs from time series
% data, based on a permutation test
%
% <matlab:open('permtest_tsdata_to_pwcgc.m') code>
%
%% Syntax
%
%    FP = permtest_tsdata_to_pwcgc(U,p,bsize,nsamps,regmode,acmaxlags,acdectol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
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
% pairwise-conditional time-domain MVGCs from the time series data |U|, based on
% randomly permuting blocks of size |bsize| of the source variable [2]. |p| is
% the model order; for other parameters see <tsdata_to_var.html |tsdata_to_var|>
% and <var_to_autocov.html |var_to_autocov|>.
%
% The first dimension of the returned matrix |FP| indexes samples, the second
% indexes the target (causee) variable and the third the source (causal)
% variable.
%
% For usage in significance testing, see <mvgc_demo_permtest.html
% |mvgc_demo_permtest|>.
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
% <permtest_tsdata_to_mvgc.html |permtest_tsdata_to_mvgc|> |
% <permtest_tsdata_to_smvgc.html |permtest_tsdata_to_smvgc|> |
% <permtest_tsdata_to_spwcgc.html |permtest_tsdata_to_spwcgc|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_pwcgc.html |autocov_to_pwcgc|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function FP = permtest_tsdata_to_pwcgc(U,p,bsize,nsamps,regmode,~,~) % v2.0 - last two parameters no longer required

if nargin < 5, regmode   = []; end % ensure default

if isempty(bsize), bsize = p; end % default to model order

[n,m,N] = size(U);
assert(m > p,'too many lags');

assert(isscalar(bsize) && isint(bsize) && bsize > 0, 'block size must be a positive integer');
nblocks  = floor(m/bsize); % number of blocks

if nblocks*bsize ~= m
    oldm = m;
    m = nblocks*bsize;
    U = U(:,1:m,:);
    fprintf(2,'WARNING: truncating sequence length by %d observation(s)\n',oldm-m);
end

FP = nan(nsamps,n,n);

for j = 1:n
    fprintf('PWCGC from node %d...',j);
    
    oj  = [1:j-1 j+1:n]; % omit j

    UPj = U;
    UUj = reshape(U(j,:,:),1,bsize,nblocks,N); % stack blocks

    for s = 1:nsamps
       %fprintf('PWCGC from node %d: permutation test sample %d of %d',j,s,nsamps); % v2.0 - reduced verbosity

        for r = 1:N
            UPj(j,:,r) = reshape(UUj(:,:,randperm(nblocks),r),1,m); % permute blocks and unstack
        end

        [A,SIG] = tsdata_to_var(UPj,p,regmode); % full regression
        if isbad(A), fprintf('\n\tERROR: VAR estimation failed'); continue; end % something went wrong
        
        [A1,C,K,res] = var_to_ss(A,SIG,2);
        if res.error, fprintf(2,'\n\t%s',res.errmsg); continue; end

        KSIGSR = K*chol(SIG,'lower');
        LSIG = log(diag(SIG));
        
        [~,SIGR,rep] = ss2iss(A1,C(oj,:),KSIGSR*KSIGSR',SIG(oj,oj),K*SIG(:,oj)); % "reduced" innovations covariance
        if rep < 0 % show-stopper!
            fprintf(2,'\n\tERROR in reduced model calculation for source node %d: ',j);
            switch rep
                case -1, fprintf(2,'DARE eigenvalues on/near unit circle');
                case -2, fprintf(2,'couldn''t find stablising DARE solution');
            end
            continue
        end
        if rep > sqrt(eps)
            fprintf(2,'\n\tWARNING in reduced model calculation for source node %d: DARE accuracy issues (relative residual = %e)\n',j,rep);
        end
        
        FP(s,oj,j) = log(diag(SIGR))-LSIG(oj); % v2.0 - loop vectorised (thanks to A. Erramuzpe Aliaga)
        
       %fprintf('\n');
    end
    fprintf('\n');
end
