%% tsdata_to_infocrit
%
% Calculate Akaike and Bayesian information criteria for VAR models from time series
% data
%
% <matlab:open('tsdata_to_infocrit.m') code>
%
%% Syntax
%
%     [aic,bic,moaic,mobic] = tsdata_to_infocrit(X,morder,regmode,verb)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
% time sries data
% _input_
%
%     X          multi-trial time series data
%     morder     maximum VAR model order or vector of model orders (see Description)
%     regmode    regression mode: 'OLS' (default) or 'LWR'
%     verb       verbosity flag: true (default) or false
%
% _output_
%
%     aic        vector of AIC values
%     bic        vector of BIC values
%     moaic      optimal model order according to AIC
%     mobic      optimal model order according to BIC
%
%% Description
%
% Calculates Akaike (AIC) and Bayesian (BIC) information criteria for VAR models
% from time series data |X|, which may be single- or multi-trial. Information
% criteria are returned in the (column) vectors |aic|, |bic|; see <infocrit.html
% |infocrit|>. Optimal model orders according to the respective criteria are
% returned in |moaic|, |mobic|.
%
% The regression mode is set by the |regmode| parameter, which may be |'LWR'|
% (default) or |'OLS'|. The former uses Morf's version of the LWR algorithm
% [1,2] while the latter calculates the OLS solution to regressions at the
% various model orders. The OLS algorithm should actually be faster, since
% regressions for all model orders are performed in "one shot" via a single QR
% decomposition; however the LWR algorithm may be slightly more accurate.
%
% The LWR algorithm must fit a VAR model at all model orders up to the maximum;
% in this case the |morder| parameter must be a scalar representing the maximum
% model order. For the OLS mode, |morder| may again be a scalar representing the
% maximum model order, or it may be a vector of model orders.
%
% If the |verb| (verbosity) flag is set (default), progress is reported;
% otherwise only warnings are reported.
%
% If VAR parameter estimation fails - because of ill-conditioned regression, or
% the residuals covariance matrix comes out as non positive-definite - a |NaN|
% is returned for that model order. This may be an indication that the time
% series data on which the autocovariance sequence estimate was based is not
% sufficiently stationary, is not of sufficient length, has (possibly lagged)
% colinearities, or has a highly skewed distribution.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] M. Morf, A. Viera, D. T. L. Lee and T. Kailath, "Recursive Multichannel
% Maximum Entropy Spectral Estimation", _IEEE Trans. Geosci. Elec._, 16(2),
% 1978.
%
%% See also
%
% <infocrit.html |infocrit|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [aic,bic,moaic,mobic] = tsdata_to_infocrit(X,morder,regmode,verb)

if nargin < 3 || isempty(regmode), regmode = 'OLS'; end % v2.0 - default changed to 'OLS', to take advantage of new algorithm
if nargin < 4 || isempty(verb),    verb    = true;  end

[n,m,N] = size(X);

assert(all(morder) < m && all(morder) > 0,'bad model order(s): must be > 0 and < %d)\n',m);

% v2.0 - in case of unforseen failure
moaic = NaN;
mobic = NaN;

X = demean(X); % no constant term

% store lags

q = max(morder);
q1 = q+1;
XX = zeros(n,q1,m+q,N);
for k = 0:q
    XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
end

if  strcmpi(regmode,'OLS') % OLS (QR decomposition) % v2.0 - new fast algorithm, based on QR decomposition

    if isscalar(morder), morder = 1:morder; end
    
    % perform the QR decomposition for all lags in one shot
    
    M  = N*(m-q);
    X0 = reshape(XX(:,1,   q1:m,:),n,  M);
    XL = reshape(XX(:,2:q1,q1:m,:),n*q,M);
    [Q,R] = qr(XL',0);
    XQ = X0*Q;

    nummo = length(morder);

    aic = nan(nummo,1);
    bic = nan(nummo,1);

    % loop through model orders

    for i = 1:nummo

        k = morder(i);

        if verb, fprintf('model order = %d',k); end
        
        nk = n*k;
        r = min(nk,M);
        wstate = warning('off','all'); lastwarn('');
        A = XQ(:,1:r)/R(1:r,1:nk)';     % OLS for regression against first k lags only
        wmsg = lastwarn; warning(wstate);        
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(A)                     % something went badly wrong
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            continue % show-stopper
        end
       
        E    = X0-A*XL(1:nk,:);        % residuals
        DSIG = det((E*E')/(M-1));      % residuals covariance matrix determinant
        if DSIG <= 0
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: residuals covariance not positive definite\n');
            continue % show-stopper
        end
        
        [aic(i),bic(i)] = infocrit(-(M/2)*log(DSIG),k*n*n,M); % no. free parameters = k*n*n, log-likelihood = -(M/2)*log(DSIG)
        if verb, fprintf(1,'\n'); end
    end

elseif strcmpi(regmode,'LWR') % LWR (Morf)

    assert(isscalar(morder),'model order must be scalar for LWR algorithm');

    q1n = q1*n;

    aic = nan(q,1);
    bic = nan(q,1);

    I = eye(n);

    % initialise recursion

    AF = zeros(n,q1n); % forward  AR coefficients
    AB = zeros(n,q1n); % backward AR coefficients (reversed compared with Morf's treatment)

    k  = 1;            % model order is k-1
    kn = k*n;
    M  = N*(m-k);
    kf = 1:kn;         % forward  indices
    kb = q1n-kn+1:q1n; % backward indices

    EF = reshape(XX(:,1:k,k+1:m,:),kn,M);
    EB = reshape(XX(:,1:k,k:m-1,:),kn,M);

    [CEF,cholp] = chol(EF*EF');
    assert(cholp == 0,'initialisation failed (''forward'' covariance matrix not positive definite)'); % v2.0 - it's a show-stopper!

    [CEB,cholp] = chol(EB*EB');
    assert(cholp == 0,'initialisation failed (''backward'' covariance matrix not positive definite)'); % v2.0 - it's a show-stopper!

    wstate = warning('off','all'); lastwarn('');
    AF(:,kf) = CEF'\I;
    AB(:,kb) = CEB'\I;
    wmsg = lastwarn; warning(wstate);        
    assert(isempty(wmsg),'initialisation failed (%s)\n',wmsg); % v2.0 - it's a show-stopper!

    % and loop

    while k <= q

        if verb, fprintf('model order = %d',k); end

        EF = AF(:,kf)*reshape(XX(:,1:k,k+1:m,:),kn,M); % forward  prediction errors
        EB = AB(:,kb)*reshape(XX(:,1:k,k:m-1,:),kn,M); % backward prediction errors

        [CEF,cholp] = chol(EF*EF');
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        [CEB,cholp] = chol(EB*EB');
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        R = CEF'\(EF*EB')/CEB;       % normalised reflection coefficients

        [CRF,cholp] = chol(I-R*R');
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        [CRB,cholp] = chol(I-R'*R);
        if cholp
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            break % it's a show-stopper!
        end

        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
        kf = 1:kn;
        kb = q1n-kn+1:q1n;

        AFPREV = AF(:,kf);
        ABPREV = AB(:,kb);

        wstate = warning('off','all'); lastwarn('');
        AF(:,kf) = CRF'\(AFPREV-R*ABPREV);
        AB(:,kb) = CRB'\(ABPREV-R'*AFPREV);
        wmsg = lastwarn; warning(wstate);        
        if ~isempty(wmsg)
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed (%s)\n',wmsg);
            break % it's a show-stopper!
        end

        wstate = warning('off','all'); lastwarn('');
        E = AF(:,1:n)\AF(:,kf)*reshape(XX(:,1:k,k+1:m,:),kn,M);
        wmsg = lastwarn; warning(wstate);        
        if ~isempty(wmsg)
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: VAR estimation failed (%s)\n',wmsg);
            break % it's a show-stopper!
        end
        
        DSIG = det((E*E')/(M-1));

        i = k-1;
        if DSIG <= 0
            if ~verb, fprintf('model order = %d',i); end
            fprintf(2,'  WARNING: residuals covariance matrix not positive definite\n');
            break % show-stopper
        end

        [aic(i),bic(i)] = infocrit(-(M/2)*log(DSIG),i*n*n,M); %  % no. free parameters = i*n*n, log-likelihood = -(M/2)*log(DSIG)i*n*n is number of parameters,  -(M/2)*log(DSIG) is log-likelihood
        if verb, fprintf(1,'\n'); end
    end

else
    error('bad regression mode ''%s''',regmode);
end

if nargout > 2 % calculate optimal model orders (NaNs are ignored)
    if isscalar(morder), morder = 1:morder; end
    [~,idx] = min(aic); moaic = morder(idx);
    [~,idx] = min(bic); mobic = morder(idx);
end
