function [LL,Lk,Lm,PACF,PACFM] = tsdata_to_varll(X,p,regmode,pacf,verb)

% NOTE: returns *AVERAGE* log-likelihood

[n,m,N] = size(X);

assert(isscalar(p) && isint(p) && p > 0 && p < m,'maximum model order must be a positive integer less than the number of observations');

if nargin < 4 || isempty(pacf), pacf = false; end
if nargin < 5 || isempty(verb), verb = false; end

p1 = p+1;
pn = p*n;
p1n = p1*n;

I = eye(n);

X = demean(X); % no constant term (don't normalise!)

% store lags

XX = zeros(n,p1,m+p,N);
for k = 0:p
    XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
end

% Note: first term is order zero!

LL = nan(p1,1); % log-likelihood
Lk = nan(p1,1); % number free parameters
Lm = nan(p1,1); % effective sample size

if pacf
	PACF = nan(n,n,p1);
	PACF(:,:,1) = I;
	PACFM = N*(m-(0:p)');
else
	PACF  = [];
	PACFM = [];
end

% order zero log-likelihood

M = N*m;

E = reshape(X,n,M); % "residuals" are just the series itself
V = (E*E')/M;       % covariance matrix (ML estimator)

LL(1) = -logdet(V)/2;
Lk(1) =  0;
Lm(1) =  M;

if  strcmpi(regmode,'OLS') % OLS

    % loop through model orders

    for k = 1:p

        if verb, fprintf('model order = %d',k); end

        k1 = k+1;
        M  = N*(m-k);

		X0 = reshape(XX(:,1,   k1:m,:),n,  M);
		XL = reshape(XX(:,2:k1,k1:m,:),n*k,M);

        wstate = warning('off','all'); lastwarn('');
        A = X0/XL;                  % OLS (QR decomposition)
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

        E = X0-A*XL;  % residuals
        V = (E*E')/M; % residuals covariance matrix (ML estimator)

        LL(k1) = -logdet(V)/2; % average LL
        Lk(k1) =  k*n*n;
        Lm(k1) =  M;

        if pacf
			A = reshape(A,n,n,k);
			PACF(:,:,k1) = A(:,:,k);
		end

        if verb, fprintf('\n'); end
    end


elseif strcmpi(regmode,'LWR') % LWR (Levinson, Wiggins & Robinson algorithm - Morf variant)

    % initialise recursion - v2.0: EF/EB Initialisation corrected (many thanks to Gonzalo Camba-Mendez for the heads-up)

    IC = inv(chol(V,'lower')); % inverse covariance square root

    k  = 1;
    kn = k*n;
    M  = N*(m-k);
    kk = 1:k;
    kf = 1:kn;         % forward  indices
    kb = p1n-kn+1:p1n; % backward indices

    AF = zeros(n,p1n); AF(:,kf) = IC; % forward  AR coefficients
    AB = zeros(n,p1n); AB(:,kb) = IC; % backward AR coefficients (reversed compared with [2])

    % LWR recursion

    while k <= p

        if verb, fprintf('model order = %d',k); end

        EF = AF(:,kf)*reshape(XX(:,kk,k+1:m,:),kn,M); % forward  prediction errors
        EB = AB(:,kb)*reshape(XX(:,kk,k:m-1,:),kn,M); % backward prediction errors

        wstate = warning('off','all'); lastwarn('');
        R = (chol(EF*EF','lower')\EF)*(chol(EB*EB','lower')\EB)'; % normalised reflection coefficients
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: reflection coefficients estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(R)                     % something went badly wrong
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: reflection coefficients estimation failed\n');
            continue % show-stopper
        end

		kp = k;
        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
		kk = 1:k;
        kf = 1:kn;
        kb = p1n-kn+1:p1n;

        AFPREV = AF(:,kf);
        ABPREV = AB(:,kb);

        wstate = warning('off','all'); lastwarn('');
        AF(:,kf) = chol(I-R*R','lower')\(AFPREV-R*ABPREV);
        AB(:,kb) = chol(I-R'*R,'lower')\(ABPREV-R'*AFPREV);
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',kp); end
            fprintf(2,'  WARNING: forward/backward VAR coefficients estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(AF) || isbad(AB) % something went badly wrong
            if ~verb, fprintf('model order = %d',kp); end
            fprintf(2,'  WARNING: forward/backward VAR coefficients estimation failed\n');
            continue % show-stopper
        end

        E = AF(:,1:n)\EF; % residuals
        V = (E*E')/M;     % residuals covariance matrix (ML estimator)

        LL(k) = -logdet(V)/2; % average LL
        Lk(k) =  kp*n*n;
        Lm(k) =  M;

        if pacf
			PACF(:,:,k) = R;
		end

        if verb, fprintf(1,'\n'); end

    end

else
    error('bad regression mode ''%s''',regmode);
end
