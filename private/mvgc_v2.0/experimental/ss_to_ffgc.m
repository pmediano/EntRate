function [F,k,dF] = ss_to_ffgc(A,C,K,V,x,y,kmax,tol)

% Infinite future (total) GC of y -> x, up to kmax steps ahead (or until convergence)

% CAREFUL! D can grow to (kmax+1)^2 x nx x n

[n,m] = ss_parms(A,C,K,V);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

if nargin < 8 || isempty(tol), tol = -Inf; end % all 'kmax' values

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
nx = length(x);
nr = length(r);
xr = 1:length(x);        % index of x in reduced quantities

L = chol(V,'lower');
V = L*L';
KL = K*L;

[KR,VR,rep] = ss2iss(A,C(r,:),KL*KL',V(r,r),K*V(:,r));
if sserror(rep), return; end
LR = chol(VR,'lower');
VR = LR*LR';
KRLR = KR*LR;

F = nan(kmax,1);

D  = L(x,:);   % D0
DR = LR(xr,:); % D0
F(1) = logdet(DR*DR')-logdet(D*D'); % F(1) is 1-step GC!
CAk1  = C(x,:);
for k = 2:kmax
	D  = [D  zeros((k-1)*nx,n ); CAk1*KL   D( end-nx+1:end,:)];
	DR = [DR zeros((k-1)*nx,nr); CAk1*KRLR DR(end-nx+1:end,:)];
	F(k) = logdet(DR*DR')-logdet(D*D');
	dF   = abs(F(k)-F(k-1));
	if dF < tol, break; end
	CAk1 = CAk1*A;
end

F = F(1:k);
