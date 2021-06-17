function [F,k,dF] = ss_to_msgc(A,C,K,V,x,y,kmax,tol)

% Finite prediction horizon GC of y -> x (k-step ahead, up to kmax or until convergence to zero)

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
xr = 1:length(x);        % index of x in reduced quantities

L = chol(V,'lower');
V = L*L';
KL = K*L;
CR = C(r,:);
[KR,VR,rep] = ss2iss(A,CR,KL*KL',V(r,r),K*V(:,r));
if sserror(rep), return; end
LR = chol(VR,'lower');
VR = LR*LR';
KRLR = KR*LR;

% CRAk1 will be CR*A^{k-1}; then BkL = CRAk1*KL is the
% the k-th MA coefficient matrix (reduced components)

F = nan(kmax,1);

Vk  = V(r,r); % k-step residuals cov mat (reduced components)
VRk = VR;     % k-step reduced residuals cov mat
F(1) = logdet(VRk(xr,xr))-logdet(Vk(xr,xr)); % F(1) is 1-step GC
CRAk1 = CR;
for k = 2:kmax
	BkL = CRAk1*KL;
	Vk = Vk + BkL*BkL';
	BRkL = CRAk1*KRLR;
	VRk = VRk + BRkL*BRkL';
	F(k)  = logdet(VRk(xr,xr))-logdet(Vk(xr,xr));
	dF    = abs(F(k));
	if dF < tol, break; end
	CRAk1  = CRAk1*A;
end

F = F(1:k);
