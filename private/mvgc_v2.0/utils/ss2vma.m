function VMAB = ss2vma(A,C,K,m)

[r,r1]  = size(A);   assert(r1 == r);
[n,r1]  = size(C);   assert(r1 == r);
[r1,n1] = size(K);   assert(n1 == n && r1 == r);

VMAB = zeros(n,n,m);
Ak = eye(size(A));
for k = 1:m
    VMAB(:,:,k) = C*Ak*K;
    Ak = A*Ak;
end
