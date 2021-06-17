function VARA = ss2var(A,C,K,m)

[r,r1]  = size(A);   assert(r1 == r);
[n,r1]  = size(C);   assert(r1 == r);
[r1,n1] = size(K);   assert(n1 == n && r1 == r);

B = A-K*C;

VARA = zeros(n,n,m);
Bk = eye(size(B));
for k = 1:m
    VARA(:,:,k) = C*Bk*K;
    Bk = B*Bk;
end
