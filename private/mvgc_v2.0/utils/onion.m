function R = onion(n,eta,Rmean)

if nargin < 2 || isempty(eta),   eta   = 1;     end % default is uniform
if nargin < 3 || isempty(Rmean), Rmean = false; end % return mean |R|

if Rmean > 0 % return theoretical mean of |R|
	nn = 1:n-1;
	f = 2*eta(:)+nn;
	R = exp(sum(nn.*(log(f-1)-log(f)),2));
	return
end

assert(nargout == 1,'too many output parameters');

R = eye(n);
b = eta+(n-2)/2;
r = 2*betarnd(b,b)-1;
R(1:2,1:2) = [1 r; r 1];
for k = 2:n-1
	b = b-1/2;
	u = rand(k,1);
	z = chol(R(1:k,1:k),'lower')*sqrt(betarnd(k/2,b))*(u/sqrt(sum(u.*u)));
	R(1:k,k+1) = z;
	R(k+1,1:k) = z';
end
