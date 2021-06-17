function [Ak,C,Kk,Vk,rep] = ss_downsample(A,C,K,V,k)

assert(isscalar(k) && isint(k) && k > 0,'Bad downsample factor (must be > 0)');

[n,r] = ss_parms(A,C,K,V);

U = K*chol(V,'lower');
KVK = U*U';
Ak = eye(r);
for m = 1:k-1
	AU = A*U;
	U = chol(AU*AU'+KVK,'lower');
    Ak = A*Ak;
end
Ak1KV = Ak*K*V; % Ak = A^(k-1)
Ak = A*Ak;      % Ak = A^k

[Kk,Vk,rep] = ss2iss(Ak,C,U*U',V,Ak1KV); % convert downsampled SS parms to ISS parms
