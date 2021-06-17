function Fms = var_to_mvgc_ms(VARA,V,x,y,q)

[n,n1,p] = size(VARA);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];              % indices of reduced variables

nx = length(x);
ny = length(y);
nr = length(r);
xr = 1:nx;              % index of x in reduced variables

pn  = p*n;
pn1 = pn-n;

Fms = nan(q,1);

% Full model ISS parameters

C = reshape(VARA,n,pn);
A = [C; eye(pn1) zeros(pn1,n)];
K = [eye(n); zeros(pn1,n)];

Vms = ISS2MSE_ms(A,C,K,V,q);      % full multi-step prediction error covariances

% Reduced model ISS parameters

AR = A;
CR = C(r,:);
KVL = K*chol(V,'lower');
[KR,VR,rep] = ss2iss(AR,CR,KVL*KVL',V(r,r),K*V(:,r));
if sserror(rep), return; end    % check DARE report, bail out on error

VRms = ISS2MSE_ms(AR,CR,KR,VR,q); % reduced multi-step prediction error covariances

% Multi-step GC

for k = 1:q
	Fms(k) = logdet(VRms(xr,xr,k)) - logdet(Vms(x,x,k));
end
