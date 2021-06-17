function H = ss2trfun(A,C,K,fres)

[n,r] = size(C);
h = fres+1;
In = eye(n);
Ir = eye(r);
H = zeros(n,n,h);
w = exp(1i*pi*((0:fres)/fres));
for k = 1:h % over [0,pi]
    H(:,:,k) = In + C*((w(k)*Ir-A)\K); % NOTE: if stable then w(k)*Ir-A always invertible!
end
