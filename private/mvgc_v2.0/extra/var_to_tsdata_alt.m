function [X,E] = var_to_tsdata_alt(A,V,W,verb)

% W must be orthogonal, zero-mean, unit-variance white noise or a scalar. In the
% latter case, Gaussian white noise of length W will be used. Note that the
% first p values are truncated, so that if a VAR of length m is required, W
% should have length p+m, where p is the number of lags.

if nargin < 4 || isempty(verb), verb = true; end

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn = p*n;
pn1 = (p-1)*n;

[n1,n2] = size(V);
assert(n1 == n2,'Residuals covariance matrix not square');
assert(n1 == n, 'Residuals covariance matrix doesn''t match VAR coefficients');

if isscalar(W)
    m = W;
    W = randn(n,m);
else
    [n1,m] = size(W);
    assert(n1 == n,'White noise doesn''t match VAR coefficients');
end
m = m-p;

X = NaN;
E = NaN;

[VSR,cholp] = chol(V,'lower');
if cholp ~= 0
    if verb, fprintf(2,'\nERROR: Residuals covariance matrix not positive-definite\n'); end
    return
end

% VAR parameters for 1-lag problem

A1 = [zeros(pn1,n) eye(pn1); reshape(flip(A,3),n,pn)];
V1 = [zeros(pn1,pn); zeros(n,pn1) V];

% Solve Lyapunov equation for 1-lag covariance matrix

try
    M1 = lyapslv('D',A1,[],-V1); % Don't balance!
catch except
    if verb, fprintf(2,'\nERROR: Failed to solve Lyapunov equation: %s\n',except.error); end
    return
end
if verb, fprintf('\nLyapunov equation relative error = %g\n',norm(M1-A1*M1*A1'-V1)/norm(M1)); end

[M1SR,cholp] = chol(M1,'lower');
if cholp ~= 0
    if verb, fprintf(2,'\nERROR: 1-lag covariance matrix not positive-definite\n'); end
    return
end

% Initial p values of X have joint covariance matrix M1

X1 = reshape(M1SR*W(1:pn)',n,p);

% Generate VAR and truncate first p values

E = VSR*W(:,p+1:p+m);
X = mvfilter([],A,[X1 E]);
X = X(:,p+1:p+m);
