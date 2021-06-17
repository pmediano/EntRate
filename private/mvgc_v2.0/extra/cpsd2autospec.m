function SS = cpsd2autospec(S)

[n,n1,h] = size(S);
assert(n1 == n,'covariance matrix not square');

SS = zeros(n,h);
for k = 1:h
	SS(:,k) = diag(S(:,:,k));
end
