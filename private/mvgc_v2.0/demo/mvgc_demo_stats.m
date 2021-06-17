%% MVGC statistics demo
%
% Demonstrates MVGC toolbox time series statistical and spectral analysis tools.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial

regmode   = [];     % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = [];     % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)
specm     = [];     % power spectrum estimation method: 'WELCH' (Welch method - default) or 'MT' (multi-taper)

etests    = false;  % do experimental (unit-root stationarity) tests? v2.0 - DEPRECATED on grounds of unreliability
stlags    = [];     % number of lags for stationarity tests (or leave empty for automatic default)

lbqmaxl   = [];     % number of lags for stationarity tests (or leave empty for automatic default)

acorr     = true;   % plot autocorrelation (else autocovariance)?

seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR data

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var5_test;
nvars = size(AT,1); % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate VAR time series data with normally distributed residuals for
% specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;

%% Unit root stationarity tests

% v2.0 DEPRECATED on grounds of unreliability
%{
if etests

    if ntrials > 1 % multitrial
        fprintf(2,'\nWARNING: unit-root stationarity tests are experimental and not really suitable for multi-trial data!\n');
    else
        fprintf(2,'\nWARNING: unit-root stationarity tests are experimental!\n');
    end

    % Augmented Dickey-Fuller unit-root test (EXPERIMENTAL)

    [adftstat,adfcval] = mvgc_adf(X,alpha,stlags);
    fprintf('\nADF statistics (critical value = %f)\n',adfcval); disp(adftstat);
    adfsig = adftstat > adfcval; % unit root; but how do we correct for multiple hypotheses?
    adfnonstat = find(adfsig);
    if isempty(adfnonstat)
        fprintf('all time series are stationary by ADF test at significance %g\n',alpha);
    else
        if ntrials > 1 % multitrial
            for r = 1:ntrials
                adfnonstat = find(adfsig(r,:));
                if ~isempty(adfnonstat)
                    fprintf(2,'WARNING: non-stationary time series by ADF test at significance %g for trial %d, variable(s): %s\n',alpha,r,num2str(adfnonstat));
                end
            end
        else
            fprintf(2,'WARNING: non-stationary time series by ADF test at significance %g for variable(s): %s\n',alpha,num2str(adfnonstat));
        end
    end

    % KPSS unit-root test (EXPERIMENTAL)

    [ksstat,kscval] = mvgc_kpss(X,alpha,stlags);
    fprintf('\nKPSS statistics (critical value = %f)\n',kscval); disp(ksstat);
    kssig = ksstat > kscval; % unit root; but how do we correct for multiple hypotheses?
    ksnonstat = find(kssig);
    if isempty(ksnonstat)
        fprintf('all time series are stationary by KPSS test at significance %g\n',alpha);
    else
        if ntrials > 1 % multitrial
            for r = 1:ntrials
                ksnonstat = find(kssig(r,:));
                if ~isempty(ksnonstat)
                    fprintf(2,'WARNING: non-stationary time series by KPSS test at significance %g for trial %d, variable(s): %s\n',alpha,r,num2str(ksnonstat));
                end
            end
        else
            fprintf(2,'WARNING: non-stationary time series by KPSS test at significance %g for variable(s): %s\n',alpha,num2str(ksnonstat));
        end
    end

end
%}

%% Model order estimation

% Calculate information criteria up to max model order

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

[~,bmo_AIC] = min(AIC);
[~,bmo_BIC] = min(BIC);

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');
xlabel('lag time');

amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
fprintf('best model order (BIC) = %d\n',bmo_BIC);
fprintf('actual model order     = %d\n',amo);

% Select model order

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = bmo_AIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = bmo_BIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation and autocovariance calculation

% Calculate VAR model; return residuals E too, since we need them later for
% statistical routines.

ptic('\n*** tsdata_to_var... ');
[A,SIG,E] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% Now calculate autocovariance according to the VAR model, to as many lags
% as it takes to decay to below the numerical tolerance level, or to acmaxlags
% lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% Report and check for errors.

var_info(info,true); % report results (and bail out on error)

% Empirical autocovariance

GE = tsdata_to_autocov(X,info.aclags);

figure(2); clf;
plot_autocov(cat(4,G,GE),{'model','data'},1/fs,[],true,acorr);
if acorr, title('Autocorrelation'); else title('Autocovariance'); end

%% Spectral analysis

ptic('*** autocov_to_cpsd... ');
[S,fres] = autocov_to_cpsd(G,fres); % for model
ptoc;

ptic('*** tsdata_to_cpsd... ');
SE = tsdata_to_cpsd(X,fres,specm);  % from data (empirical)
ptoc;

% plot (auto-)spectra

figure(3); clf;
plot_cpsd(cat(4,S,SE),{'model','data'},fs,[],true);
title('Autospectra');

%% VAR stats tests

% Check that residuals are white (Durbin-Watson test).

[dw,dwpval] = whiteness(X,E);
fprintf('\nDurbin-Watson statistics:\n'); disp(dw);
dwsig = significance(dwpval,alpha,mhtc); % significance adjusted for multiple hypotheses
notwhite = find(dwsig);
if isempty(notwhite)
    fprintf('all residuals are white by Durbin-Watson test at significance %g\n',alpha);
else
    fprintf(2,'WARNING: autocorrelated residuals at significance %g for variable(s): %s\n',alpha,num2str(notwhite));
end

% Check that residuals are white (Ljung-Box portmanteau test).

if isempty(lbqmaxl), lbqmaxl = ceil(log(nobs)); end

Q = lbqtest(E,morder,lbqmaxl);
Qc = lbqtest_cval(nvars,morder,lbqmaxl,alpha);

figure(4); clf;
plot_tsdata([Q Qc]',{'Q','Q(crit)'},1/fs);
title('Ljung-Box portmanteau test');
legend('Q','Qc','location','southeast');
xlabel('lag time');

fprintf('\nLjung-Box statistics:\n');
fprintf('lag : %s\n',num2str(morder+1:lbqmaxl,       '%9d'));
fprintf('Q   : %s\n',num2str(Q(morder+1:lbqmaxl)',   '%9.2f'));
fprintf('Qc  : %s\n\n',num2str(Qc(morder+1:lbqmaxl)','%9.2f'));

if Q(lbqmaxl) < Qc(lbqmaxl);
    fprintf('residuals are white by Ljung-Box portmanteau test (%d lags), at significance %g\n',lbqmaxl,alpha);
else
    fprintf(2,'WARNING: residuals NOT white by Ljung-Box portmanteau test (%d lags) at significance %g\n',lbqmaxl,alpha);
end

% Check R^2 stats.

[~,RSQADJ] = rsquared(X,E);
fprintf('\nRSQ (adjusted):\n'); disp(RSQADJ);
rsqthreshold = 0.3; % like GCCA
badqsq = find(RSQADJ < rsqthreshold);
if isempty(badqsq)
    fprintf('adjusted r-squares OK: > %g%% of variance is accounted for by the model\n',100*rsqthreshold);
else
    fprintf(2,'WARNING: low adjusted r-square values (< %g) for variable(s): %s\n',rsqthreshold,num2str(badqsq));
end

% Check model consistency (ie. proportion of correlation structure of the data
% accounted for by the model).

cons = 100*consistency(X,E); % percent
fprintf('\nmodel consistency = %.0f%%\n',cons);
consthreshold = 80;          % like GCCA
if cons > consthreshold
    fprintf('consistency OK: > %g%%\n',consthreshold);
else
    fprintf(2,'WARNING: low consistency (< %g%%)\n',consthreshold);
end

%%
% <mvgc_demo_stats.html back to top>
