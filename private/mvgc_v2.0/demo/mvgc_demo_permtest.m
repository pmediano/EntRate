%% MVGC permutation test demo
%
% Demonstrates permutation significance testing with MVGC on generated VAR data
% for a 5-node network with known causal structure (see <var5_test.html
% |var5_test|>). Pairwise-conditional Granger causalities are estimated and
% significance tested using both the theoretical and permutation test null
% distributions. The theoretical and permutation test null distributions for all
% pairs are plotted together and compared using a Kolmogorov-Smirnov test (see
% <matlab:doc('kstest') |kstest|>).
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
% _Aust. N. Z. J. Stat._, 43(1), 2001.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 10;     % number of trials
nobs      = 100;    % number of observations per trial
nperms    = 100;    % number of permutations for permutation test
bsize     = [];     % permutation test block size: empty for automatic (uses model order)

regmode   = [];     % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)

alpha     = 0.05;   % significance level for all statistical tests
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR data

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var5_test;
nvars = size(AT,1);  % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate VAR time series data with normally distributed residuals for
% specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;

%% VAR model estimation and autocovariance calculation

morder = size(AT,3); % actual model order; on real data - i.e. with no generative model
                     % available - use information criteria to estimate (see 'mvgc_demo')

% Calculate VAR model

ptic('*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

%% Convert VAR to SS

% Convert the VAR model to innovations-form SS model.

fprintf('*** var_to_ss\n');
[A1,C,K,ssinfo] = var_to_ss(A,SIG); % A1 is the "companion matrix" of A

assert(~ssinfo.error,'Bailing out'); % abort on error

%% Granger causality estimation

% Calculate time-domain pairwise-conditional causalities from SS model.

ptic('*** ss_to_pwcgc... ');
F = ss_to_pwcgc(A1,C,K,SIG);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Theoretical significance test (adjusting for multiple hypotheses).

pval_t = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2);
sig_t  = significance(pval_t,alpha,mhtc);

%% Permutation test

ptic('\n*** tsdata_to_mvgc_pwc_permtest\n');
FNULL = permtest_tsdata_to_pwcgc(X,morder,bsize,nperms);
ptoc('*** tsdata_to_mvgc_pwc_permtest took ',[],1);

% (We should probably check for failed permutation estimates here.)

% Permutation test significance test (adjusting for multiple hypotheses).

pval_p = empirical_pval(F,FNULL);
sig_p  = significance(pval_p,alpha,mhtc);

%% Plot causalities, p-values, significance and cdfs

figure(1); clf;
subplot(2,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(2,3,2);
plot_pw(pval_t);
title({'p-values';'(theoretical)'});
subplot(2,3,3);
plot_pw(sig_t);
title({['Significant at p = ' num2str(alpha)];'(theoretical)'});
subplot(2,3,5);
plot_pw(pval_p);
title({'p-values';'(perm test)'});
subplot(2,3,6);
plot_pw(sig_p);
title({['Significant at p = ' num2str(alpha)];'(perm test)'});

% Plot empirical vs theoretical null cdfs, perform Kolmogorov-Smirnov tests on
% null/theoretical distributions.

figure(2); clf;
pval_ks = nan(nvars); % KS test p-value
fmax = max(FNULL(:));
k = 0;
for i = 1:nvars
    for j = 1:nvars
        k = k+1;
        if i ~= j
            FN = sort(FNULL(:,i,j));
            PN = empirical_cdf(FN,FNULL(:,i,j));
            PT = 1-mvgc_pval(FN,morder,nobs,ntrials,1,1,nvars-2); % theoretical cdf = 1 - p-value
            subplot(nvars,nvars,k);
            plot(FN,[PN PT]);
            title(sprintf('%d -> %d',j,i));
            xlim([0 fmax]);
            [~, pval_ks(i,j)] = kstest(FN,[FN PT]); % 2-sided KS test of empirical vs. theoretical null distribution
        end
    end
end

% KS significance test (adjusting for multiple hypotheses).

sig_ks = significance(pval_ks,alpha,mhtc); % 1 = distributions are significantly different at level alpha

fprintf('\nKS-test p-values:\n\n'); disp(pval_ks);
fprintf('KS-test significance: ones indicates that permutation and theoretical\nnull distributions are significantly different at significance level %g; \n\n',alpha); disp(sig_ks);

%%
% <mvgc_demo_permtest.html back to top>
