function B = VAR2VMA(VARA,q)

[n,n1,p] = size(VARA);
assert(n1 == n,'VAR coefficients matrix has bad shape');

pn = p*n;
pn1 = pn-n;

C = reshape(VARA,n,pn);
A = [C; eye(pn1) zeros(pn1,n)]; % the 1-lag "companion matrix"
K = [eye(n); zeros(pn1,n)];

B = zeros(n,n,q);
Ak = eye(pn);
for k = 1:q
	B(:,:,k) = C*Ak*K;
	Ak = A*Ak;
end
