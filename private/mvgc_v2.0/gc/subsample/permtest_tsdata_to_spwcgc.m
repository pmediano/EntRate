%% permtest_tsdata_to_spwcgc
%
% Calculate null distribution for pairwise-conditional frequency-domain MVGCs
% from time series data, based on a permutation test
%
% <matlab:open('permtest_tsdata_to_spwcgc.m') code>
%
%% Syntax
%
%    fP = permtest_tsdata_to_spwcgc(U,p,fres,bsize,nsamps,regmode,acmaxlags,acdectol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     p          model order (number of lags)
%     fres       frequency resolution (default: automatic)
%     bsize      permutation block size (default: use model order)
%     nsamps     number of permutations
%     regmode    regression mode (default as for 'tsdata_to_var')
%     acmaxlags  maximum autocovariance lags  (default as for 'var_to_autocov')
%     acdectol   autocovariance decay tolerance (default as for 'var_to_autocov')
%
% _output_
%
%     fP         permutation test spectral Granger causalities (null distribution)
%
%% Description
%
% Returns |nsamps| samples from the empirical null distribution of the
% pairwise-conditional frequency-domain MVGCs from the time series data |U|,
% based on randomly permuting blocks of size |bsize| of the source variable [2].
% |p| is the model order; for other parameters see <tsdata_to_var.html
% |tsdata_to_var|> and <var_to_autocov.html |var_to_autocov|>.
%
% The first dimension of the returned matrix |fP| indexes samples, the second
% indexes the target (causee) variable, the third the source (causal)
% variable and the fourth frequency.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
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
% <permtest_tsdata_to_pwcgc.html |permtest_tsdata_to_pwcgc|> |
% <permtest_tsdata_to_smvgc.html |permtest_tsdata_to_smvgc|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_spwcgc.html |autocov_to_spwcgc|> |
% <sfreqs.html |sfreqs|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function fP = permtest_tsdata_to_spwcgc(U,p,fres,bsize,nsamps,regmode,~,~) % v2.0 - last two parameters no longer required

% v2.0 - completely rewritten using a new algorithm; doesn't need
% 'autocov_xform()' anymore. See autocov_to_smvgc.m.

if nargin < 6, regmode   = []; end % ensure default

if isempty(bsize), bsize = p; end % default to model order

[n,m,N] = size(U);
assert(m > p,'too many lags');

assert(isscalar(bsize) && isint(bsize) && bsize > 0, 'block size must be a positive integer');
nblocks  = floor(m/bsize); % number of blocks

if nblocks*bsize ~= m
    oldm = m;
    m = nblocks*bsize;
    U = U(:,1:m,:);
    fprintf(2,'WARNING: truncating sequence length by %d observations\n',oldm-m);
end

h = fres+1;

fP = nan(nsamps,n,n,h);

for j = 1:n
    fprintf('spectral PWCGC from node %d...',j);
    
    oj  = [1:j-1 j+1:n]; % omit j

    UPj = U;
    UUj = reshape(U(j,:,:),1,bsize,nblocks,N); % stack blocks

   for s = 1:nsamps
       %fprintf('spectral PWCGC from node %d: permutation test sample %d of %d',j,s,nsamps); % v2.0 - reduced verbosity

        for r = 1:N
            UPj(j,:,r) = reshape(UUj(:,:,randperm(nblocks),r),1,m); % permute blocks and unstack
        end

        [A,SIG] = tsdata_to_var(UPj,p,regmode); % full regression
        if isbad(A), fprintf('\n\tERROR: VAR estimation failed'); continue; end % something went wrong
        
        [A1,C,K,res] = var_to_ss(A,SIG,2);
        if res.error, fprintf(2,'\n\t%s',res.errmsg); continue; end
        
        KSIGSR = K*chol(SIG,'lower');

        CR = C(oj,:);
        [KR,SIGR,rep] = ss2iss(A1,CR,KSIGSR*KSIGSR',SIG(oj,oj),K*SIG(:,oj)); % "reduced" innovations covariance
        if rep < 0 % show-stopper!
            fprintf(2,'\n\tERROR in reduced model calculation for source node %d: ',j);
            switch rep
                case -1, fprintf(2,'DARE eigenvalues on/near unit circle\n');
                case -2, fprintf(2,'couldn''t find stablising DARE solution\n');
            end
            continue
        end
        if rep > sqrt(eps)
            fprintf(2,'\n\tWARNING in reduced model calculation for source node %d: DARE accuracy issues (relative residual = %e)\n',j,rep);
        end
        
        H = ss2trfun(A1,C,K,fres);

        BR = ss2itrfun(A1,CR,KR,fres); % reduced inverse transfer function

        for ii = 1:n-1;
            i  = oj(ii);         % i index in omitted j indices
            oi = [1:i-1 i+1:n];  % omit i

            SR  = SIGR(ii,ii);   % reduced spectrum is flat!
            LSR = log(SR);

            PSIGSR = chol(parcov(SIG,oi,i),'lower'); % partial covariance square root

            for k = 1:h
                HRk = BR(ii,:,k)*H(oj,oi,k)*PSIGSR; % transformed reduced transfer function
                fP(s,i,j,k) = LSR - log(SR-HRk*HRk');
            end
        end

       %fprintf('\n');
    end
    fprintf('\n');
end
