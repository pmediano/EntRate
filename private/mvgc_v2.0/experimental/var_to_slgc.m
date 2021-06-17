function [F,info] = var_to_slgc(A,V,x,y,p,verb)

% Single-lag GC of y -> x at 1:p lags

if nargin < 6 || isempty(verb), verb = false; end

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

[G,info] = var_to_autocov(A,V);

q1 = size(G,3);
q = q1-1;
qn = q*n;
nx = length(x);
u = 1:n;

if p > q % pad G with zeros
	G = cat(3,G,zeros(n,n,p-q));
end

% Construct full autocovariance matrix

GG = zeros(qn,qn);               % full autocovariance matrix
for i = 0:q-1
	for j = i:q-1
		GG(u+i*n,u+j*n) = G(:,:,j-i+1);
	end
end
GG = symmetrise(GG);
[~,cholp] = chol(GG,'lower');
assert(cholp == 0,'Full autocovariance matrix not positive-definite');


C = G(x,x);
g = reshape(G(x,n+1:end),nx,qn); % G1 ... Gq
LDV = logdet(V(x,x));

F  = nan(p,1);
for k = 1:p
	if verb, fprintf('lag = %2d of %2d\n',k,p); end

	o = 1:qn; o((k-1)*n+y) = []; % omit just k-th lag of y

	L = g(:,o)/chol(GG(o,o));

	F(k) = logdet(C-L*L')-LDV;
end
