function [A,C,K,V,Z,E] = tsdata_to_ss(X,pf,r)

% Estimate an (innovations form) state space model from an empirical observation
% time series using Larimore's Canonical Correlations Analysis (CCA) state
% space-subspace algorithm (W. E. Larimore, in Proc. Amer. Control Conference,
% vol. 2, IEEE, 1983).
%
% X         - observation process time series
% pf        - past/future horizons for canonical correlations
% r         - SS model order: 'SVC' (default), 'prompt', or a positive integer
%
% A,C,K,V   - estimated innovations form state space parameters
% Z         - estimated state process time series
% E         - estimated innovations process time series
%
% The past/future horizons pf may be supplied as a 2-vector [p,f] or a scalar p
% = f = pf. Bauer recommends setting p = f = 2*p, where p is the optimal
% VAR model order for the observation process X according to Aikaike's Information
% Criterion (AIC).

[n,m,N] = size(X);

assert(all(isint(pf(:))),'past/future horizon must be a 2-vector or a scalar positive integer');
if isscalar(pf)
    p = pf;    f = pf;
elseif isvector(pf) && length(pf) == 2
    p = pf(1); f = pf(2);
else
    error('past/future horizon must be a 2-vector or a scalar positive integer');
end

A = NaN;
C = NaN;
K = NaN;
V = NaN;
Z = NaN;
E = NaN;

X = demean(X); % no constant term (don't normalise!)

mp  = m-p;
mp1 = mp+1;
mf  = m-f;
mpf = mp1-f; % m-p-f+1

M  = N*mp;
M1 = N*mp1;
Mh = N*mpf;

Xf = zeros(n,f,mpf,N);
for k = 1:f
    Xf(:,k,:,:) = X(:,p+k:mf+k,:);
end
Xf = reshape(Xf,n*f,Mh);

XP = zeros(n,p,mp1,N);
for k = 0:p-1
    XP(:,k+1,:,:) = X(:,p-k:m-k,:);
end
Xp = reshape(XP(:,:,1:mpf,:),n*p,Mh);
XP = reshape(XP,n*p,M1);

[Wf,cholp] = chol((Xf*Xf')/Mh,'lower');
assert(cholp == 0,'forward weight matrix not positive definite');

[Wp,cholp] = chol((Xp*Xp')/Mh,'lower');
assert(cholp == 0,'backward weight matrix not positive definite');

BETA = Xf/Xp; % 'OH' estimate: regress future on past
assert(all(isfinite(BETA(:))),'subspace regression failed');

[~,S,U] = svd(Wf\BETA*Wp); % SVD of CCA-weighted OH estimate
assert(all(isfinite(S(:))),'SVD failed');

sval = diag(S); % the singular values

Z = reshape((diag(sqrt(sval(1:r)))*U(:,1:r)'/Wp)*XP,r,mp1,N); % Kalman states estimate; note that Z starts at t = p+1, has length mp1 = m-p+1
assert(all(isfinite(Z(:))),'Kalman states estimation failed');

% Calculate model parameters by regression

C = reshape(X(:,p+1:m,:),n,M)/reshape(Z(:,1:mp,:),r,M);       % observation matrix
assert(all(isfinite(C(:))),'C parameter estimation failed');

E = reshape(X(:,p+1:m,:),n,M) - C*reshape(Z(:,1:mp,:),r,M);   % innovations
V = (E*E')/(M-1);                                             % innovations covariance matrix (unbiased estimate)

AK = reshape(Z(:,2:mp1,:),r,M)/[reshape(Z(:,1:mp,:),r,M);E];
assert(all(isfinite(AK(:))),'A/K parameter estimation failed');

A = AK(:,1:r);                                                % state transition matrix
K = AK(:,r+1:r+n);                                            % Kalman gain matrix

E = reshape(E,n,mp,N);
