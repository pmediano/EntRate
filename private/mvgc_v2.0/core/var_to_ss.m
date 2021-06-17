% READY TO GO!

function [A,C,K,info] = var_to_ss(VARA,V,report)

[n,n1,p] = size(VARA);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn1 = (p-1)*n;

C = reshape(VARA,n,p*n);
A = [C; eye(pn1) zeros(pn1,n)]; % the 1-lag "companion matrix"
K = [eye(n); zeros(pn1,n)];

if nargout > 3
    assert(nargin > 1,'Need covariance matrix for VAR info');
    if nargin < 3, report = []; end % use ss_info default
    info = ss_info(A,C,K,V,report);
end
