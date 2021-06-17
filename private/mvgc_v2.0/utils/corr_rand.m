function [R,L,retries,iters] = corr_rand(n,g,vexp,iscorr,tol,maxretries,maxiters)

% Generate random correlation matrix with specified generalised correlation (or
% "multi-information" - see below). Creates a uniform random orthogonal matrix
% and a random variance vector with independent chi^2(1) distribution, then
% adjusts the variance vector until within tolerance of specified generalised
% correlation.
%
% Generalised correlation is supplied as negative log-generalised correlation or
% "multi-information" g = -log|R|, where R is the correlation matrix and n the
% number of dimensions. For n = 2, g = -log(1-rho^2) where rho is the Pearson
% correlation coefficient (so for small rho, g ~ rho^2).
%
% If g is not supplied (or is empty) a correlation matrix is sampled uniformly
% at random over the space of n x n correlation matrices using the "onion"
% algorithm. If g is zero, the identity matrix is returned
%
% If the 'iscorr' flag is set, then g is considered a "correlation":
% rho = sqrt(1-exp(-(2/n)*g)) and is converted to g = -(n/2)*log(1-rho^2).
%
% NOTE: if you want a covariance matrix with standard deviations in the (column)
% vector s, then use V = s.*R.*s', or L -> s.*L
%
% n          - number of dimensions
% g          - multi-information (g = -log|R|): g = 0 yields zero correlation
% vexp       - variance exponent (see 'corr_rand_exponent')
% iscorr     - g is supplied as "correlation" rather than multi-information
% tol        - numerical tolerance (default: sqrt(eps))
% maxretries - maximum retries to find large enough correlation (default: 1000)
% maxiters   - maximum iterations for binary chop (default: 1000)
%
% R          - correlation matrix
% L          - Cholesky (left) factor: L*L' = R
% retries    - number of retries required
% iters      - number of binary chop iterations required

if nargin < 2                         g          = [];        end
if nargin < 3 || isempty(vexp),       vexp       = 2;         end
if nargin < 4 || isempty(iscorr),     iscorr     = false;     end
if nargin < 5 || isempty(tol),        tol        = sqrt(eps); end
if nargin < 6 || isempty(maxretries), maxretries = 1000;      end
if nargin < 7 || isempty(maxiters),   maxiters   = 1000;      end

L       = NaN;
retries = 0;
iters   = 0;

if isempty(g)
    R = onion(n,-g); % random on space of n x n correlation matrices (uniform if g = -1; see 'onion.m')
    [L,pchol] = chol(R,'lower');
    if pchol ~= 0
        fprintf(2,'ERROR: ''corr_rand'' result not positive-definite (onion failed)\n');
        R = NaN;
        L = NaN;
        return
    end
    R = L*L';
    return
end

assert(isscalar(g) && isnumeric(g),'(log-generalised) correlation must be empty or a non-negative scalar');

if g < eps && g > -eps
	R = eye(n);
	L = eye(n);
    return
end

assert(g > 0,'(log-generalised) correlation must be empty or a non-negative scalar');

if iscorr % g is supplied as "correlation" r = sqrt(1-exp((2/n)*g))
	g = -(n/2)*log(1-g*g);
end

% We calcluate a (positive-definite) covariance matrix with given generalised
% correlation, and finally convert it to a correlation matrix. If V is a
% variance-covariance matrix, then
%
% g = sum(log(diag(V)))-logdet(V))

gtarget = g; % target value for g

% Find random orthogonal M and vector v of variances such that for V = M*diag(v)*M'
% (which will be pos-def), g is >= gtarget. The rationale is that adding the same
% constant c to all variances always decreases g, so we may use a binary chop to
% home in on gtarget. Note that since M is orthogonal, |V| = prod(v), so that we
% may calculate g efficiently as sum(log(diag(V))-log(v)).

gotit = false;
for retries = 0:maxretries
    [Q,R] = qr(randn(n));
	v = realpow(abs(randn(n,1)),vexp);
    M = Q*diag(sign(diag(R))); % M orthogonal
    V = M*diag(v)*M';          % V is pos-def
    g = sum(log(diag(V))-log(v));
    if g >= gtarget
        gotit = true;
        break
    end
end
if ~gotit
    fprintf(2,'ERROR: ''corr_rand'' timed out on retries (g too large?)\n');
    R = NaN;
    return
end
D = diag(V);

% Got M and v (for the binary chop we just need the variances D of V and the
% unrotated variances v). Now set binary chop initial high value so that g <
% gtarget; start c at 1 and keep doubling

c = 1;
while g > gtarget
    g = sum(log(D+c) - log(v+c));
    c = 2*c;
end
chi = c;
clo = 0;
ghi = gtarget+tol;
glo = gtarget-tol;

% Do binary chop

gotit = false;
for iters = 1:maxiters % binary chop
    c = (clo+chi)/2;
    g = sum(log(D+c) - log(v+c));
    if     g < glo % g too small, set chi to current
        chi = c;
    elseif g > ghi % g too big, set clo to current
        clo = c;
    else
        gotit = true;
        break;
    end
end
if ~gotit % this shouldn't really happen :-/
    fprintf(2,'ERROR: ''corr_rand'' timed out on binary-chop (numerical instability?)\n');
    R = NaN;
    return
end

% We have a c that meets tolerance

V = M*diag(v+c)*M';

% Check V really is pos-def

[L,pchol] = chol(V,'lower');
if pchol ~= 0
    fprintf(2,'ERROR: ''corr_rand'' result not positive-definite\n');
    R = NaN;
    L = NaN;
    return
end

% Convert to correlation matrix

L = diag(1./sqrt(sum(L.*L,2)))*L;
R = L*L';
