function B = mvconv(A,v)

[n,m,p] = size(A);
p = p-1;
q = length(v)-1;
r = p+q;

B = zeros(n,m,r+1);
for k = 0:r
	Bk = zeros(n,m);
	for i = max(0,k-q):min(p,k)
		Bk = Bk + A(:,:,i+1)*v(k-i+1);
	end
	B(:,:,k+1) = Bk;
end

%{
% can use conv()

A = permute(A,[3,1,2]);

BB = zeros(n,m,r+1);
for i = 1:n
	for j = 1:n
		BB(i,j,:) = conv(A(:,i,j),v);
	end
end

maxabs(B-BB)
%}
