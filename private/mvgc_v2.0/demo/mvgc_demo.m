%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test data generation

ntrials   = 4;       % number of trials
nobs      = 500;     % number of observations per trial
fs        = 200;     % sample rate (Hz)

% Actual VAR model generation parameters

tnet      = tnet5;   % connectivity network
moact     = 6;       % model order
rho       = 0.95;    % spectral radius
wvar      = 0.5;     % var coefficients decay weighting factor
rmi       = 0.5;     % residuals log-generalised correlation (multi-information)
                     % g = -log|R|. g = 0 yields zero correlation,g = [] is uniform random
                     % on space of correlation matrices

% VAR model order estimation

moregmode = 'LWR';   % VAR model estimation regression mode ('OLS' or 'LWR')
mosel     = 'LRT';   % model order selection ('ACT', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
momax     = 2*moact; % maximum model order for model order selection

% VAR model parameter estimation

regmode   = 'LWR';   % VAR model estimation regression mode ('OLS' or 'LWR')

% MVGC (time domain) statistical inference

tstats    = 'dual';  % test statistic ('single', 'dual' or 'both')
alpha     = 0.05;    % significance level for Granger casuality significance test
mhtc      = 'FDR';   % multiple hypothesis test correction (see routine 'significance')

% MVGC (frequency domain)

fres      = [];      % spectral MVGC frequency resolution (empty for automatic calculation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('seed', 'var'), seed  = 0; end % random seed (0 for unseeded)
if ~exist('plotm','var'), plotm = 0; end % plot mode (figure number offset, or Gnuplot terminal string)

%% Generate random VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

% Seed random number generator.

rng_seed(seed);

% Generate random VAR coefficients for test network.

AA = var_rand(tnet,moact,rho,wvar);
nvars = size(AA,1); % number of variables

% Generate random residuals covariance (in fact correlation) matrix.

VV = corr_rand(nvars,rmi);

% Report information on the generated VAR model and check for errors.

infoo = var_info(AA,VV);
assert(~infoo.error,'VAR error(s) found - bailing out');

% Generate multi-trial VAR time series data with normally distributed residuals
% for generated VAR coefficients and residuals covariance matrix.

ptic('*** varfima_to_tsdata... ');
X = varfima_to_tsdata(AA,[],[],VV,nobs,ntrials);
ptoc;

% Remove temporal mean and normalise by temporal variance.
% Not strictly necessary, but may help numerical stability
% if data has very large or very small values.

X = demean(X,true);

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate and plot VAR model order estimation criteria up to specified maximum model order.

ptic('\n*** tsdata_to_varmo... ');
if isnumeric(plotm), plotm = plotm+1; end
[moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,momax,moregmode,[],[],plotm);
ptoc;

% Select and report VAR model order.

morder = moselect(sprintf('VAR model order selection (max = %d)',momax), mosel,'ACT',moact,'AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
assert(morder > 0,'selected zero model order! GCs will all be zero!');
if morder >= momax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,V] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed - bailing out');

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A,V);
assert(~info.error,'VAR error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Estimated time-domain pairwise-conditional Granger causalities

ptic('*** var_to_pwcgc... ');
[F,stats] = var_to_pwcgc(A,V,tstats,X,regmode);
ptoc;
assert(~isbad(F,false),'GC estimation failed');

% Significance test (F- and likelihood ratio), adjusting for multiple hypotheses.

sigF  = significance(stats.(tstats).F.pval, alpha,mhtc);
sigLR = significance(stats.(tstats).LR.pval,alpha,mhtc);

% For comparison, we also calculate the actual pairwise-conditional causalities

ptic('*** var_to_pwcgc... ');
FF = var_to_pwcgc(AA,VV);
ptoc;
assert(~isbad(FF,false),'GC calculation failed');

% Plot time-domain causal graph, p-values and significance.

maxF = 1.1*max(nanmax(F(:),nanmax(FF(:))));
pdata = {FF,F;sigF,sigLR};
ptitle = {'PWCGC (actual)','PWCGC (estimated)'; sprintf('F-test (%s-regression)',tstats),sprintf('LR test (%s-regression)',tstats)};
maxp = [maxF maxF;1 1];
if isnumeric(plotm), plotm = plotm+1; end
plot_gc(pdata,ptitle,[],maxp,plotm);

%% Granger causality estimation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% Calculate spectral pairwise-conditional causalities resolution from VAR model
% parameters. If not specified, we set the frequency resolution to something
% sensible (based on the spectral radii of the VAR model - see var_info) - we also
% warn if the calculated resolution is very large, as this may cause problems.

if isempty(fres)
    fres = 2^nextpow2(max(info.acdec,infoo.acdec)); % alternatively, fres = 2^nextpow2(nobs);
	fprintf('\nUsing frequency resolution %d\n',fres);
end
if fres > 10000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

ptic(sprintf('\n*** var_to_spwcgc (at frequency resolution = %d)... ',fres));
f = var_to_spwcgc(A,V,fres);
ptoc;
assert(~isbad(f,false),'spectral GC estimation failed');

% For comparison, we also calculate the actual pairwise-conditional spectral causalities

ptic(sprintf('*** var_to_spwcgc (at frequency resolution = %d)... ',fres));
ff = var_to_spwcgc(AA,VV,fres);
ptoc;
assert(~isbad(ff,false),'spectral GC calculation failed');

% Get frequency vector according to the sampling rate.

freqs = sfreqs(fres,fs);

% Plot spectral causal graphs.

if isnumeric(plotm), plotm = plotm+1; end
plot_sgc({ff,f},freqs,'Spectral Granger causalities (blue = actual, red = estimated)',plotm);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities. Note that this may occasionally fail if a certain condition
% on the VAR parameters is not satisfied (see refs. [4,5]).

Fint = bandlimit(f,3); % integrate spectral MVGCs (frequency is dimension 3 of CPSD array

fprintf('\n*** GC spectral integral check... ');
rr = abs(F-Fint)./(1+abs(F)+abs(Fint)); % relative residuals
mrr = max(rr(:));                       % maximum relative residual
if mrr < 1e-5
    fprintf('PASS: max relative residual = %.2e\n',mrr);
else
    fprintf(2,'FAIL: max relative residual = %.2e (too big!)\n',mrr);
end

%%
% <mvgc_demo.html back to top>
