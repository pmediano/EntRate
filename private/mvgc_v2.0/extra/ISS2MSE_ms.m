function W = ISS2MSE_ms(A,C,K,V,q)

[r,r1]  = size(A); assert(r1 == r);
[n,r1]  = size(C); assert(r1 == r);
[r1,n1] = size(K); assert(n1 == n && r1 == r);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

L = K*chol(V,'lower');
W = zeros(n,n,q);

CAk1 = C;
W(:,:,1) = V;
for k = 1:q-1
	CAk1L = CAk1*L;
	W(:,:,k+1) = W(:,:,k) + CAk1L*CAk1L';
	CAk1 = CAk1*A;
end
