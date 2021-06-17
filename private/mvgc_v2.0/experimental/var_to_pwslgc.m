function [F,info] = var_to_pwslgc(A,V,p,verb)

% Single-lag GC of y -> x at 1:p lags

if nargin < 4 || isempty(verb), verb = false; end

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

[G,info] = var_to_autocov(A,V);

q1 = size(G,3);
q = q1-1;
qn = q*n;
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

C = diag(G(u,u));
g = reshape(G(u,n+1:end),n,qn); % G1 ... Gq
LV = log(diag(V));

F = nan(p,n,n);
for y = 1:n
	if verb, fprintf('source = %2d of %2d ',y,n); end
	r = [1:y-1,y+1:n];
	for k = 1:p
		if verb, fprintf('.'); end

		o = 1:qn; o((k-1)*n+y) = []; % omit just k-th lag of y

		L = g(r,o)/chol(GG(o,o));

		F(k,r,y) = log(C(r)-diag(L*L'))-LV(r);
	end
end
