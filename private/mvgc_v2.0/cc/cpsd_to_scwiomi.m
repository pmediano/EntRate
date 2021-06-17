function C = cpsd_to_scwiomi(S)

[n,n1,h] = size(S);
assert(n1 == n,'CPSD must be square at each frequency');

LDS = zeros(h,1);
for k = 1:h
    LDS(k)  = logdet(S(:,:,k));
end

C = nan(n,h);
for i = 1:n
	oi = 1:n; oi(i) = [];
    for k = 1:h
		C(i,k) = log(S(i,i,k)) + logdet(S(oi,oi,k)) - LDS(k);
	end
end
