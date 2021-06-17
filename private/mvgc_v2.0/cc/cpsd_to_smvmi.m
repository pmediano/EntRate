function C = cpsd_to_smvmi(S,x,y)

[n,n1,h] = size(S);
assert(n1 == n,'CPSD must be square at each frequency');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z = 1:n; z([x y]) = [];
xz = [x,z];
yz = [y,z];

C = nan(1,h);

for k = 1:h
	C(k) = logdet(S(xz,xz,k)) + logdet(S(yz,yz,k)) - logdet(S(z,z,k)) - logdet(S(:,:,k));
end
