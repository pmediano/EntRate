function A = VEC2VAR(PI,G)

[n,~,p1] = size(G);
p = p1+1;

A = zeros(n,n,p);
A(:,:,1) = PI + eye(n) + G(:,:,1);
for k = 2:p1
    A(:,:,k) = G(:,:,k) - G(:,:,k-1);
end
A(:,:,p) = -G(:,:,p1);
